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
from XICRA.scripts import pilfer_caller
from XICRA.config import set_config
from XICRA.modules import database

import HCGB.functions.fasta_functions as HCGB_fasta
import HCGB.functions.time_functions as HCGB_time
import HCGB.functions.system_call_functions as HCGB_sys
import HCGB.functions.files_functions as HCGB_files
import HCGB.functions.main_functions as HCGB_main
import HCGB.functions.aesthetics_functions as HCGB_aes
import HCGB.format_conversion.file_splitter as HCGB_splitter

##########################################################3
def annotate_sam_call(sam_file, gold_piRNA, ncpu, folder, Debug):
    
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

    if Debug:
        print("seq_dict gold piRNA")
        ##print(seq_dict)
        print ("ATTENTION: Very big file. See file provided for details: " + gold_piRNA)

    ## as it might be very big, we are splitting and processing in parallel
    path_given = HCGB_files.create_folder(os.path.join(folder, "split_sam"))
    HCGB_splitter.split_file_call(sam_file, ncpu*8, "split_file", False, 'SAM', os.path.abspath(path_given), Debug)
    
    list_sam_files = HCGB_main.get_fullpath_list(path_given, Debug)
    
    ## remove non-desired files
    list_sam_files = [s for s in list_sam_files if '.parsed' not in s]
    list_sam_files = [s for s in list_sam_files if '.split_file_success' not in s]

    ## send for each subset using multiple threads
    cpus2use=ncpu*4 ## Each single process uses just 10% each cpu, we would increase speed by using x4 times more. 
		    ## It might not be best way but it really speeds up the process and does not consum many RAM or CPUs
    with concurrent.futures.ThreadPoolExecutor(max_workers=cpus2use) as executor:
        commandsSent = { executor.submit(pilfer_caller.annotate_sam, seq_id, 
                                         subset_sam, Debug): subset_sam for subset_sam in list_sam_files }

        for cmd2 in concurrent.futures.as_completed(commandsSent):
            details = commandsSent[cmd2]
            try:
                data = cmd2.result()
            except Exception as exc:
                print ('***ERROR:')
                print (cmd2)
                print('%r generated an exception: %s' % (details, exc))

    list_sam_files = (s + '.parsed' for s in list_sam_files)
    
    if Debug:
        print(list_sam_files)

    return(list_sam_files)

################################
def process_call(bam_file, sample_folder, name, gold_piRNA, ncpu, Debug):
    
    ## check if previously trimmed and succeeded
    filename_stamp = sample_folder + '/.success_BAM2PILFER'
    if os.path.isfile(filename_stamp):
        print ("\n+ Converting BAM file into PILFER input file")
        stamp = HCGB_time.read_time_stamp(filename_stamp)
        print (colored("\tA previous command generated results on: %s [%s -- %s]" %(stamp, name, 'BAMtoPILFER'), 'yellow'))
        
        pilfer_file = os.path.join(sample_folder, name + ".pilfer.bed")
        return (pilfer_file)

    else:
        code_returned = bam2pilfer(bam_file, sample_folder, name, gold_piRNA, ncpu, Debug)
        if code_returned:
            HCGB_time.print_time_stamp(filename_stamp)
            return (code_returned) ## file name
        else:
            print ('** Sample %s failed...' %name)
            return(False)

################################
def merge_sam_bed(bed_file, sam_file, pilfer_tmp, Debug):

    ## TODO:
    ## Linux only: find an alternative
    ## paste Aligned.sortedByCoord.out.bed Aligned.sortedByCoord.out.sam | awk -v "OFS=\t" '{print $1, $2, $3, $16, $6}'

    if Debug:
        HCGB_aes.debug_message("bed_file: " + bed_file, "yellow")
        HCGB_aes.debug_message("sam_file: " + sam_file, "yellow")
        HCGB_aes.debug_message("pilfer_tmp: " + pilfer_tmp, "yellow")

    cmd_paste = "paste %s %s | grep \':P\' | awk -v \"OFS=\t\" \'{print $1, $2, $3, $16, $6}\' > %s" %(bed_file, sam_file, pilfer_tmp)
        
    paste_code = HCGB_sys.system_call(cmd_paste, False, True)
    return paste_code

################################
def bam2pilfer(bam_file, out_folder, name, annot_info, ncpu, Debug):
    
    out_folder = os.path.abspath(out_folder)
    
    print("+ Convert BAM to PILFER input file")
    print("+ Several pre-processing steps are required:")
    
    ## convert BAM2bed
    print("\t- Convert original BAM to BED format file...")
    bed_file = bedtools_caller.convert_bam2bed(name + '_original', bam_file, out_folder, pilfer=True, debug=Debug)
    
    ## substract ncRNA included (no piRNA)
    print("\t- Subtract ncRNA annotated reads from BAM to reduce processing...")
    bed_file_reduced = bedtools_caller.subtract_coordinates(bed_file, annot_info['general']['ncRNA'], 
                                                            out_folder, name + '_subtracted', "", Debug)
    
    ## get IDs and convert original bam to sam only with retained IDs
    ## required samtools 1.12
    print("\t- Get subtracted reads...")
    ## generate paste filter tmp file
    reads_to_get = os.path.join(out_folder, "reads2retain.txt")
    filename_stamp = out_folder + '/.subtract_read_success'
    if os.path.isfile(filename_stamp) and HCGB_files.is_non_zero_file(reads_to_get):
        stamp = HCGB_time.read_time_stamp(filename_stamp)
        print (colored("\tA previous command generated results on: %s [%s -- %s]" %(stamp, name, 'subtract reads'), 'yellow'))
    else:
        
        original_stdout = sys.stdout # Save a reference to the original standard output
        with open(reads_to_get, 'w') as f_out:
            sys.stdout = f_out # Change the standard output to the file we created.
            for line in open(bed_file_reduced):
                line=line.rstrip().split("\t")
                read=line[3].split("/")[0]
                print(read)

        sys.stdout = original_stdout
        f_out.close()
    
    ## convert bamtosam
    print("\t- Convert BAM to SAM format file only including subtracted reads...")
    sam_file = samtools_caller.bam2sam(name + '_reduced', bam_file, out_folder, ncpu, "-h -N " + reads_to_get, False, Debug)
    #                                sample, input_file, path_given, ncpu, options, sam2bam=False, Debug=False):

    ## convert sam to bam
    print("\t- Convert SAM to BAM format file only including subtracted reads...")
    bam_file_reduced = samtools_caller.bam2sam(name + '_reduced', sam_file, out_folder, ncpu, "", True, Debug)

    print("\t- Convert reduced BAM to BED format file...")
    bed_file_mapping = bedtools_caller.convert_bam2bed(name + '_reduced', bam_file_reduced, out_folder, pilfer=True, debug=Debug)
    
    ## New files
    pilfer_tmp = os.path.join(out_folder, name + ".pilfer.bed.tmp")
    pilfer_file = os.path.join(out_folder, name + ".pilfer.bed")
    
    ## Annotate reads using gold piRNA
    print("\t- Annotate reads in BAM using gold piRNA information provided...")
    
    ## generate paste filter tmp file
    filename_stamp = out_folder + '/.annotate_sam_success'
    if os.path.isfile(filename_stamp):
        stamp = HCGB_time.read_time_stamp(filename_stamp)
        print (colored("\tA previous command generated results on: %s [%s -- %s]" %(stamp, name, 'annotate sam'), 'yellow'))
        
        path_given = HCGB_files.create_folder(os.path.join(os.path.abspath(out_folder), "split_sam"))
        list_sam_parsed = HCGB_main.get_fullpath_list(path_given, Debug)
   
        ## remove non-desired files
        list_sam_parsed = [s for s in list_sam_parsed if not '.parsed' not in s]
        list_sam_parsed = [s for s in list_sam_parsed if '.split_file_success' not in s]

    else:
        list_sam_parsed = annotate_sam_call(sam_file, annot_info['piRBase']['gold_piRNA'], ncpu, out_folder, Debug)
        ## print time stamp
        HCGB_time.print_time_stamp(filename_stamp)
    
    ## merge annotated sam
    parsed_sam = os.path.join(out_folder, name + '_reduced-annotated.sam')
    list_sam_parsed = sorted(list_sam_parsed, key=lambda x: int("".join([i for i in x if i.isdigit()])))

    if Debug:
        HCGB_aes.debug_message("List of sam files parsed:")
        print(list_sam_parsed)

    ## add annotated reads into asingle file
    parsed_sam_file = open(parsed_sam, "w")
    for l in list_sam_parsed:
        for row in open(l, 'r'):
            parsed_sam_file.write(row)

    parsed_sam_file.close()
    
    #####
    print("\t- Create PILFER format file...")
    filename_stamp = out_folder + '/.convert_bam2pilfer_success'
    if os.path.isfile(filename_stamp):
        if (HCGB_files.is_non_zero_file(pilfer_tmp)):
            print (colored("\tA previous command generated results on: %s [%s -- %s]" %(stamp, name, 'bam2pilfer'), 'yellow'))
            return(pilfer_file)

    ## generate paste filter tmp file
    code_exec = merge_sam_bed(bed_file_mapping, parsed_sam, pilfer_tmp, Debug)
    if not code_exec:
        print ("** Some error ocurred while merging annotated SAM and BED file")
        exit()

    ## create bed file summarized:
    ## cat Aligned.sortedByCoord.out.sam.pilfer.bed.tmp | bedtools groupby -g 1,2,3,4,5 -c 4 -o count 

    bedtools_exe = set_config.get_exe("bedtools", Debug)
    cmd_bedtools = "%s sort -chrThenSizeA -i %s | sort -k 4 | %s groupby -o count -g 1,2,3,4,5 -c 4 | awk -v \"OFS=\t\" \'{print $1, $2, $3, $4, $6, $5}\' > %s " %(bedtools_exe, pilfer_tmp, bedtools_exe,  pilfer_file)
    bed_code = HCGB_sys.system_call(cmd_bedtools, False, True)
    if not bed_code:
        print("** Some error occurred while generating PILFER input file")
        exit()

    ## remove tmp files
    #os.remove(pilfer_tmp)
    #os.remove(sam_file_out)

    HCGB_time.print_time_stamp(filename_stamp)
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

    parser.add_argument('--database', 
                        help='Absolute path for database containing piRNA subfolder.', required=True);    

    parser.add_argument("-t", "--threads", type=int, help="Number of CPUs to use [Default: 2].", default=2)

    parser.add_argument("-s", "--species", help="Species abbreviation. [Default: hsa].", default="hsa")


    args=parser.parse_args();

    
    args.path_given = HCGB_files.create_folder(args.path_given)
    
    args.database = os.path.abspath(args.database)
    args.input = os.path.abspath(args.input)
    
    ## retrieved piRNA information
    annot_info = database.piRNA_info(args.database, args.species, False)
    
    ### 
    process_call(args.input, args.path_given, args.name, annot_info, args.threads, False)
    

################################
if __name__== "__main__":
    main()
