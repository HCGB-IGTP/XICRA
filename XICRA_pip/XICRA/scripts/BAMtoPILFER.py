#usr/bin/env python

## useful imports
import time
import io
import os
import re
import sys
from sys import argv
import subprocess

from XICRA.scripts import bedtools_caller

################################
def process_call(bam_file, sample_folder, name, gold_piRNA, Debug):
    
    ## check if previously trimmed and succeeded
    filename_stamp = sample_folder + '/.success'
    if os.path.isfile(filename_stamp):
        stamp = functions.time_functions.read_time_stamp(filename_stamp)
        print (colored("\tA previous command generated results on: %s [%s -- %s]" %(stamp, name, 'BAMtoPILFER'), 'yellow'))
    else:
        print ("\n+ Converting BAM file into PILFER input file")

        code_returned = conversion(bam_file, out_folder, name, gold_piRNA, Debug)
        if code_returned:
            functions.time_functions.print_time_stamp(filename_stamp)
        else:
            print ('** Sample %s failed...' %name)
            
##########################################################
def get_known_piRNA(fasta_file, Debug):
    """
    
    """
    
    fasta_count = defaultdict(int)
    with open(fasta_file, 'r') as fh:
        lines = []
        for line in fh:
            lines.append(line.rstrip())
            if len(lines) == 2:
                record = HCGB_fasta.process_fasta(lines)
                print(record)
                #sys.stderr.write("Record: %s\n" % (str(record)))
                lines = []
                fasta_count[record['sequence']] = record['name'].split(">")[1]
                
    return fasta_count        

################################
def conversion(bam_file, out_folder, name, gold_piRNA_file, Debug):
    
    ## get known piRNA seqs
    gold_piRNA = get_known_piRNA(gold_piRNA_file, Debug)
    
    ## convert BAM2bed
    bed_file = bedtools_caller.convert_bam2bed(name, bam_file, out_folder, pilfer=True, debug=Debug)
    
    ## convert bamtosam
    sam_file = samtools_caller.conversion(bam_file, "SAM", options="no_header", Debug)

    ## Variables
    dirname_name = os.path.dirname(bam_file)
    split_name = os.path.splitext( os.path.basename(bam_file) )
    bed_file = folder + '/' + split_name[0] + '.bed'
    sam_file = folder + '/' + split_name[0] +  '.sam'
    pilfer_tmp = folder + '/' + split_name[0] + '.tmp.pilfer.bed'
    pilfer_file = folder + '/' + split_name[0] + '.pilfer.bed'
    
    
    ## generate samtools
    if (os.path.isfile(sam_file)):
        print ("\t+ File %s already exists" %sam_file)
    else:
        cmd_samtools = "%s view %s > %s" %(samtools_exe, bam_file, sam_file)
        output_file.write(cmd_samtools)
        output_file.write("\n")        
        try:
            subprocess.check_output(cmd_samtools, shell = True)
        except Exception as exc:
            print ('***ERROR:')
            print (cmd_samtools)
            print('samtools view command generated an exception: %s' %exc)
            exit()
    
    ## generate paste filter tmp file
    if (os.path.isfile(pilfer_tmp)):
        print ("\t+ File %s already exists" %pilfer_tmp)
    else:
        ## paste Aligned.sortedByCoord.out.bed Aligned.sortedByCoord.out.sam | awk -v "OFS=\t" '{print $1, $2, $3, $16, $6}' 
        cmd_paste = "paste %s %s | awk -v \"OFS=\t\" \'{print $1, $2, $3, $16, $6}\' > %s" %(bed_file, sam_file, pilfer_tmp)
        output_file.write(cmd_paste)
        output_file.write("\n")        
        try:
            subprocess.check_output(cmd_paste, shell = True)
        except Exception as exc:
            print ('***ERROR:')
            print (cmd_paste)
            print('paste bed sam command generated an exception: %s' %exc)
            exit()
    
    ## parse pilfer tmp file
    counter = 1
    previous_line = ()
    
    # Open file OUT
    output_file = open(pilfer_file, 'w')
    
    # Open file IN
    fileHandler = open (pilfer_tmp, "r")
    while True:
        # Get next line from file
        line = fileHandler.readline().strip()
           # If line is empty then end of file reached
        if not line :
            break;
    
        seq = line.split('\t')[3]
        real_seq = seq.split('::PU')
        seq_len = len(str(real_seq[0]))
    
        ## Discard smaller
        if (previous_line):
            if (previous_line == line):
                line = previous_line
                counter += 1
            else:
                line_split = previous_line.split('\t')
                output_file.write('%s\t%s\t%s\t%s::PI\t%s\t%s\n' %(line_split[0], line_split[1], line_split[2], line_split[3], counter, line_split[4]))
    
        #counter += 1
        while True:
            #get next line
            next_line = fileHandler.readline().strip()
        
            if (next_line == line):
                counter += 1
            else:
                line_split = line.split('\t')
                output_file.write('%s\t%s\t%s\t%s::PI\t%s\t%s\n' %(line_split[0], line_split[1], line_split[2], line_split[3], counter, line_split[4]))
    
                previous_line = next_line
                counter = 1
                break;
 
    ## close and finish
    fileHandler.close()
    output_file.close()

################################
def main():
    ## this code runs when call as a single script
    parser=argparse.ArgumentParser(description='''Create piRNA analysis using PILFER''');
    
    parser.add_argument('--input', '-i', help='Input BAM file', required=True);

    parser.add_argument('--name', '-n',
                        help='Name of the sample. Default: use filename provided.', default="");

    parser.add_argument('--known_piRNA', 
                        help='Absolute path for piRBase gold piRNA sequences.', required=True);    

################################
if __name__== "__main__":
    main()
