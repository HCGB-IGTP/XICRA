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

##########################################################
def pilfer_module_call(sample_folder, name, bam_file, database_folder, threads, species, Debug):
    print()
    
    ## check if previously trimmed and succeeded
    filename_stamp = sample_folder + '/.success'
    if os.path.isfile(filename_stamp):
        stamp = functions.time_functions.read_time_stamp(filename_stamp)
        print (colored("\tA previous command generated results on: %s [%s -- %s]" %(stamp, name, 'pilfer'), 'yellow'))
    else:
        
        ## retrieved piRNA information
        annot_info = database.piRNA_info(database_folder, species, Debug)
        code_returned = pilfer_caller(sample_folder, name, bam_file, annot_info, threads, Debug)
        if code_returned:
            functions.time_functions.print_time_stamp(filename_stamp)
        else:
            print ('** Sample %s failed...' %name)

##########################################################
def pilfer_caller(sample_folder, name, bam_file, annot_info, threads, Debug):
    """
    
    """
    
    ## convert BAM to PILFER Input file
    bam_pilfer = BAMtoPILFER.process_call(bam_file, sample_folder, name, annot_info.gold_piRNA, threads, True)

    ## create clusters using pilfer.py
    
    
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