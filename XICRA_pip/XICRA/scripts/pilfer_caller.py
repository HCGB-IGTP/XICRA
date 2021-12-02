#!/usr/bin/env python3
############################################################
## Jose F. Sanchez                                        ##
## Copyright (C) 2019-2021 Lauro Sumoy Lab, IGTP, Spain   ##
############################################################
## useful imports
import time
import io
import os
import re
import sys
from sys import argv
import subprocess
import argparse
from termcolor import colored
import pandas as pd
from collections import defaultdict

import numpy as np
from operator import itemgetter
import sys, csv, os, argparse

from XICRA.config import set_config
from XICRA.scripts import bedtools_caller
from XICRA.scripts import BAMtoPILFER
from XICRA.modules import database

import HCGB.functions.aesthetics_functions as HCGB_aes
import HCGB.functions.files_functions as HCGB_files
import HCGB.functions.system_call_functions as HCGB_sys
import HCGB.functions.fasta_functions as HCGB_fasta
import HCGB.functions.time_functions as HCGB_time
import HCGB.format_conversion

##########################################################3
def pilfer_code(infile, outfile):
    """
    Predicts piRNA clusters from an input BED file
    
    This code is copy from PILFER software. See copyright and License details in https://github.com/rishavray/PILFER
    Original code: June 2018
    https://github.com/rishavray/PILFER/blob/master/tools/pilfer.py
    
    Modifications: November 2021
    - Add some comments to understand code and clarify it.
    - Add to read BED file 
    
    """
    
    ## Input example
    ##1    997178    997194    CGGGTTCATTTCCCGGATGATGCACCA::PU    1    +
    ##1    1312232    1312248    GGCGAATAGTCGGTGGCTGTAGACGCGAAACCT::PU    1    +
    ##1    1355449    1355465    ATTCGCGGATGAGCAGCGAAATGCGCAGCGTC::PU    1    +
    ##1    1909803    1909819    TTGAATTTGCCGAACCGGCTTTATGGAACC::PU    2    +
    ##1    2527004    2527020    GTGCAGATCGGCACGCCAGGGCGAAT::PU    1    -
    ##1    3473028    3473044    GCCCGACCACCAGGCCAGACCGTCAAT::PU    1    -
    ##1    3707624    3707640    AGCACGCCGCCAGCGTTCATCCTGAG::PU    5    -
    ##1    3772285    3772301    TCGGCCCACGGGCTCGTGGCGTCGGA::PU    3    +
    ##1    3772450    3772466    ATTCCGGACGAGGCGTCGCTCGATCAA::PU    1    +
    ##1    4353091    4353107    TTTCGCTTGAGGACGCTACCGATTTCATTC::PU    1    +
    
    #Variables
    sd_factor=3 ## "The factor by which the read should be away from standard deviation to be called a peak"
    count_reads = []
    chrom_dict = {}
    results = []

    infile_pt = open(infile, "rb")
    #csvin = csv.reader(infile_pt, delimiter = "\t")
    
    #Reading the BED records
    for row in infile_pt:
        field=row.strip().split('\t')
        count_reads.append(float(row[4]))
        
        ## classify in a dictionary by reference sequence
        if field[0] in chrom_dict:
            chrom_dict[field[0]].append(field)
        else:
            chrom_dict[field[0]] = [field]
        
    #Mean and SD calculations
    mean = np.mean(count_reads)
    sd = np.std(count_reads)
    
    #Sorting the dictionary
    for key in chrom_dict:
        chrom_dict[key] = sorted(chrom_dict[key],key=itemgetter(1))
    #del chrom_dict['chrY']
    
    #print "Chromosome\tStart\tEnd\tScore"
    #Cluster calculation
    for key in chrom_dict:
        j = 0
        while j < len(chrom_dict[key]):
            row = chrom_dict[key][j]
            if (row[4] - mean)/sd >= sd_factor:
                start_bp = row[2] - 100000
                if start_bp < 0:
                    start_bp = 0
                end_bp = row[2]
                score = 0
                new_score = 0
                start_read_index = -1
                end_read_index = 0
                
                cur_read_index = chrom_dict[key].index(row)
    
                #Calculate the initial score and start index
                for read in chrom_dict[key]:
                    if read[1] >= start_bp and read[1] <= start_bp + 100000:
                        score += read[4]
                        if start_read_index == -1:
                            start_read_index = chrom_dict[key].index(read)
                            start_bp = read[1]
                        end_read_index = chrom_dict[key].index(read)
    
                new_score = score
                max_start = start_read_index
                max_end = end_read_index
                #calculating optimum 100KB window
                ##for i in xrange(start_read_index+1,cur_read_index+1): ## xrange not available in python3
                for i in range(start_read_index+1,cur_read_index+1): 
                    new_score = score - chrom_dict[key][i-1][4]
                    while  end_read_index+1 < len(chrom_dict[key]) and chrom_dict[key][end_read_index +1][1] <= chrom_dict[key][i][1] + 100000 :
                        new_score += chrom_dict[key][end_read_index+1][4]
                        end_read_index += 1
                    
                    if new_score > score:
                        score = new_score
                        max_start = i
                        max_end = end_read_index
    
                #print (key + ":" + str(chrom_dict[key][max_start][1]) + "-" + str(chrom_dict[key][max_end][2]) + "\t" + str(score)) ## pilfer default format
                string2print = key + "\t" + str(chrom_dict[key][max_start][1]) + "\t" + str(chrom_dict[key][max_end][2]) + "\t" + str(score) ## bed format
                results.append(string2print)
                j = max_end
            j += 1

    ## print output file
    outfile_pt = open(outfile, "w")
    for l in results:
        outfile_pt.write(l)
    outfile_pt.close()
    
    ## return information in bed format
    return (results)
            
##########################################################3
def annotate_sam(seq_id, sam_file, Debug):
    """
    This code is copy from PILFER software. See copyright and License details in https://github.com/rishavray/PILFER
    Original code: June 2018
    https://github.com/rishavray/PILFER/blob/master/tools/annotate_sam.py
    
    Modifications: November 2021
    - Add some comments to understand code and clarify it.
    - Add to read fasta file as dict and not as list as provided. 
    
    This particular scripts uses input in sam format and a file with gold piRNA sequences to identify putative and well-knonw piRNAs.
    
    """
    
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

        ## take into account sam header
        if (field[0].startswith('@')):
            #sam_file_write.write(row)
            continue

        seq = field[9]
        
        if (int(field[1]) & (0x10)):
            seq = HCGB_fasta.ReverseComplement(seq)
        
        ## Set putative (PU), known piRNA (PI) or none
        if seq in seq_id:
            field.append("XP:Z:PI")
            field[9] = field[9] + '::PI'
        #elif len(row[9])>=mini and len(field[9])<=maxi and int(f[0])>=100:
        elif len(field[9])>=mini and len(field[9])<=maxi:
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

##########################################################
def pilfer_module_call(sample_folder, name, bam_file, database_folder, threads, species, Debug):
    print()
    
    ## check if previously trimmed and succeeded
    filename_stamp = sample_folder + '/.success'
    if os.path.isfile(filename_stamp):
        stamp = HCGB_time.read_time_stamp(filename_stamp)
        print (colored("\tA previous command generated results on: %s [%s -- %s]" %(stamp, name, 'pilfer'), 'yellow'))
    else:
        
        ## retrieved piRNA information
        annot_info = database.piRNA_info(database_folder, species, Debug)
        code_returned = pilfer_caller(sample_folder, name, bam_file, annot_info, threads, Debug)
        if code_returned:
            HCGB_time.print_time_stamp(filename_stamp)
        else:
            print ('** Sample %s failed...' %name)

##########################################################
def pilfer_caller(sample_folder, name, bam_file, annot_info, threads, Debug):
    """
    Create required input files, calls pilfer to create clusters and intersect with transposon information coordinates.
    """
    
    ## convert BAM to PILFER Input file
    bam_pilfer = BAMtoPILFER.process_call(bam_file, sample_folder, name, annot_info, threads, Debug)
    
    ## create clusters using pilfer.py
    outfile = os.path.join(sample_folder, "pilfer_clusters.bed")
    ## Call pilfer
    pilfer_code(bam_pilfer, outfile)
    
    ## intersect with Transposon data
    bedtools_caller.intersect_coordinates(bam_pilfer, annot_info["TEsmall_db"]["TE"], sample_folder, 'transposon_intersect', '', Debug)
    
    
    return()


#######################################################
def main():
    ## this code runs when call as a single script
    parser=argparse.ArgumentParser(description='''Create piRNA analysis using PILFER''');
    
    parser.add_argument('--input', '-i', help='Input BAM file', required=True);

    parser.add_argument('--name', '-n',
                        help='Name of the sample. Default: use filename provided.', default="");
    
    parser.add_argument('--path_given', '-p',
                        help='Path to save results generated. Default: use path from filename provided.', default="");

    parser.add_argument('--database', '-db',
                        help='Absolute path for XICRA database containing piRNA files.', default="");

    parser.add_argument("-t", "--threads", type=int, help="Number of CPUs to use [Default: 2].", default=2)

    args=parser.parse_args();
    
    ## Lets convert bam2bed
    args.path_given = HCGB_files.create_folder(args.path_given)
    
    ## lets split the big file provided
    files_generated = pilfer_module_call(args.path_given, args.name, args.input, 
                                         args.database, args.threads, "hsa", Debug=True)
    
    return ()


######
if __name__== "__main__":
    main()
