#usr/bin/env python

## useful imports
import time
import io
import os
import re
import sys
import csv
from sys import argv
import subprocess
from collections import defaultdict
import argparse
from termcolor import colored
import concurrent.futures

from XICRA.scripts import bedtools_caller
from XICRA.scripts import samtools_caller
from XICRA.config import set_config

import HCGB.functions.fasta_functions as HCGB_fasta
import HCGB.functions.time_functions as HCGB_time
import HCGB.functions.system_call_functions as HCGB_sys
import HCGB.functions.files_functions as HCGB_files
import HCGB.functions.main_functions as HCGB_main

import HCGB.format_conversion.file_splitter as HCGB_splitter

##########################################################3
def annotate_sam(seq_id, sam_file, Debug):
           
    ## create output
    sam_file_out = sam_file + ".parsed"
    sam_file_write = open(sam_file_out, "w")
    
    fileReader = open(sam_file, 'r')
    
    ## Give a minimum and maximum length
    mini = 26
    maxi = 33
    
    #csvin = csv.reader(fileReader, delimiter="\t")
    #csvout = csv.writer(sam_file_write, delimiter="\t")
    for row in open(sam_file):
        ## Original
        #f = row[0].split(":")
        #row[0] = f[1]
        field=row.strip().split('\t')
        seq = field[9]
        
        if (int(field[1]) & (0x10)):
            seq = HCGB_fasta.ReverseComplement(seq)
        
        ## Set putative (PU), known piRNA (PI) or none
        if seq in seq_id:
            field.append("XP:Z:PI")
            field[9] = field[9] + '::PI'
        #elif len(row[9])>=mini and len(field[9])<=maxi and int(f[0])>=100:
        elif len(row[9])>=mini and len(field[9])<=maxi:
            field.append("XP:Z:PU")
            field[9] = field[9] + '::PU'
        #else:
        ## too big or not known piRNA 
        
        ## append length in all
        field.append("XC:i:"+str(len(seq)))
    
        sam_file_write.write("\t".join(field) + "\n")
        
    ## close files
    sam_file_write.close()
    fileReader.close()

##########################################################3
def annotate_sam_call(sam_file, sam_file_out, gold_piRNA, ncpu, folder, Debug):
    
    """
    This code is copy from PILFER software. See copyright and License details in https://github.com/rishavray/PILFER
    Original code: June 2018
    https://github.com/rishavray/PILFER/blob/master/tools/annotate_sam.py
    
    Modifications: November 2021
    - Add some comments to understand code and clarify it.
    - Add to read fasta file as dict and not as list as provided. 
    
    This particular scripts uses input in sam format and a file with gold piRNA sequences to identify putative and well-knonw piRNAs.
    
    """

    ### Original code
    # def ReverseComplement(seq):
    #    seq_dict = {'A':'T','T':'A','G':'C','C':'G','N':'N'}
    #    return "".join([seq_dict[base] for base in reversed(seq)])
    #
    # seq_id = []
    # with open(sys.argv[1],"rb") as listin:
    #    for seq in listin:
    #        seq_id.append(seq.strip())
    # seq_id = list(set(seq_id))
    
    
    ## read fasta file and save results in list
    seq_dict = HCGB_fasta.get_fasta_dict(gold_piRNA, Debug=Debug)
    seq_id = list(seq_dict.keys())
    
    ## as it might be very big, we are splitting and processing in parallel
    path_given = HCGB_files.create_folder(os.path.join(folder, "split_sam"))
    HCGB_splitter.split_file_call(sam_file, ncpu*10, "spli_file", False, 'SAM', os.path.abspath(path_given), Debug)
    
    list_sam_files = HCGB_main.get_fullpath_list(path_given, Debug)
    
    ## send for each subset using multiple threads
    with concurrent.futures.ThreadPoolExecutor(max_workers=ncpu*2) as executor:
        commandsSent = { executor.submit(annotate_sam, seq_id, subset_sam, Debug): subset_sam for subset_sam in list_sam_files }

        for cmd2 in concurrent.futures.as_completed(commandsSent):
            details = commandsSent[cmd2]
            try:
                data = cmd2.result()
            except Exception as exc:
                print ('***ERROR:')
                print (cmd2)
                print('%r generated an exception: %s' % (details, exc))

    return(True)

################################
def process_call(bam_file, sample_folder, name, gold_piRNA, ncpu, Debug):
    
    ## check if previously trimmed and succeeded
    filename_stamp = sample_folder + '/.success'
    if os.path.isfile(filename_stamp):
        print ("\n+ Converting BAM file into PILFER input file")
        stamp = HCGB_time.read_time_stamp(filename_stamp)
        print (colored("\tA previous command generated results on: %s [%s -- %s]" %(stamp, name, 'BAMtoPILFER'), 'yellow'))
    else:
        code_returned = bam2pilfer(bam_file, sample_folder, name, gold_piRNA, ncpu, Debug)
        if code_returned:
            HCGB_time.print_time_stamp(filename_stamp)
        else:
            print ('** Sample %s failed...' %name)


################################
def bam2pilfer(bam_file, out_folder, name, gold_piRNA_file, ncpu, Debug):
    
    
    print("+ Convert BAM to PILFER input file")
    print("+ Several pro-processing steps are required:")
    
    ## convert BAM2bed
    print("\t- Convert BAM to BED format file...")
    bed_file = bedtools_caller.convert_bam2bed(name, bam_file, out_folder, pilfer=True, debug=Debug)
    
    ## convert bamtosam
    print("\t- Convert BAM to SAM format file...")
    sam_file = samtools_caller.bam2sam(name, bam_file, out_folder, ncpu, header=False, Debug=Debug)
    sam_file_out = sam_file + ".tmp"
    
    ## New files
    pilfer_tmp = os.path.join(out_folder, os.path.basename(sam_file) + ".pilfer.bed.tmp")
    pilfer_file = os.path.join(out_folder, os.path.basename(sam_file) + ".pilfer.bed")
    
    ## Annotate reads using gold piRNA
    print("\t- Annotate reads in BAM using gold piRNA information provided...")
    
    ## generate paste filter tmp file
    filename_stamp = out_folder + '/.annotate_sam_success'
    if HCGB_files.is_non_zero_file(sam_file_out) and os.path.isfile(filename_stamp):
        stamp = HCGB_time.read_time_stamp(filename_stamp)
        print (colored("\tA previous command generated results on: %s [%s -- %s]" %(stamp, name, 'annotate sam'), 'yellow'))
    else:
        annotate_sam_call(sam_file, sam_file_out, gold_piRNA_file, ncpu, out_folder, Debug)
        ## print time stamp
        HCGB_time.print_time_stamp(filename_stamp)
    
    ## generate paste filter tmp file
    ## paste Aligned.sortedByCoord.out.bed Aligned.sortedByCoord.out.sam | awk -v "OFS=\t" '{print $1, $2, $3, $16, $6}'
    
    print("\t- Create PILFER format file...")
    
    filename_stamp = out_folder + '/.convert_bam2pilfer_success'
    if os.path.isfile(filename_stamp):
        if (HCGB_files.is_non_zero_file(pilfer_tmp)):
            print (colored("\tA previous command generated results on: %s [%s -- %s]" %(stamp, name, 'bam2pilfer'), 'yellow'))
            return(pilfer_file)

    ## TODO:
    ## Linux only: find an alternative
    
    cmd_paste = "paste %s %s | awk -v \"OFS=\t\" \'{print $1, $2, $3, $16, $6}\' > %s" %(bed_file, sam_file_out, pilfer_tmp)
    
    paste_code = HCGB_sys.system_call(cmd_paste, False, True)
    if paste_code:
        ## print time stamp
        HCGB_time.print_time_stamp(filename_stamp)
    else:
        print("Some error ocurred during the conversion from BAM to PILFER input...")
        exit()
        
    ## parse pilfer tmp file
    
    ## create bed file summarized:
    ## cat Aligned.sortedByCoord.out.sam.pilfer.bed.tmp | bedtools groupby -g 1,2,3,4,5 -c 4 -o count 

    bedtools_exe = set_config.get_exe("bedtools", Debug)
    cmd_bedtools = "cat %s | %s groupby -o count -g 1,2,3,4,5 -c 4 | grep '::P'> %s " %(pilfer_tmp, bedtools_exe,  pilfer_file)
    bed_code = HCGB_sys.system_call(cmd_bedtools, False, True)
    if not bed_code:
        print("** Some error occurred while generating PILFER input file")
        exit()

    ## remove tmp files
    os.remove(pilfer_tmp)
    os.remove(sam_file_out)
 
    return(pilfer_file)    
 
################################
def main():
    ## this code runs when call as a single script
    parser=argparse.ArgumentParser(description='''Create piRNA analysis using PILFER''');
    
    parser.add_argument('--input', '-i', help='Input BAM file', required=True);

    parser.add_argument('--name', '-n',
                        help='Name of the sample. Default: use filename provided.', default="");

    parser.add_argument('--path_given', '-p',
                        help='Name of the path to store results. Default: use filename provided.', default="");

    parser.add_argument('--gold_piRNA', 
                        help='Absolute path for piRBase gold piRNA sequences.', required=True);    

    parser.add_argument("-t", "--threads", type=int, help="Number of CPUs to use [Default: 2].", default=2)


    args=parser.parse_args();

    
    args.path_given = HCGB_files.create_folder(args.path_given)
    
    ### 
    process_call(args.input, args.path_given, args.name, args.gold_piRNA, args.threads, True)
    

################################
if __name__== "__main__":
    main()
