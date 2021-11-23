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

import HCGB.functions.aesthetics_functions as HCGB_aes
import HCGB.functions.files_functions as HCGB_files
import HCGB.functions.system_call_functions as HCGB_sys
import HCGB.functions.fasta_functions as HCGB_fasta
import HCGB.functions.time_functions as HCGB_time
import HCGB.format_conversion

##########################################################
def pilfer_module_call(sample_folder, name, bed_sample, annot_bed, Debug):
    print()
    
    ## check if previously trimmed and succeeded
    filename_stamp = sample_folder + '/.success'
    if os.path.isfile(filename_stamp):
        stamp = functions.time_functions.read_time_stamp(filename_stamp)
        print (colored("\tA previous command generated results on: %s [%s -- %s]" %(stamp, name, 'pilfer'), 'yellow'))
    else:
        code_returned = pilfer_caller(sample_folder, name, bed_sample, annot_bed, Debug)
        if code_returned:
            functions.time_functions.print_time_stamp(filename_stamp)
        else:
            print ('** Sample %s failed...' %name)

##########################################################
def pilfer_caller(sample_folder, name, bed_sample, annot_bed, Debug):
    print()


#######################################################
def main():
    ## this code runs when call as a single script
    parser=argparse.ArgumentParser(description='''Create piRNA analysis using PILFER''');
    
    parser.add_argument('--input', '-i', help='Input BAM file', required=True);

    parser.add_argument('--name', '-n',
                        help='Name of the sample. Default: use filename provided.', default="");
    
    parser.add_argument('--path_given', '-p',
                        help='Path to save results generated. Default: use path from filename provided.', default="");

    parser.add_argument('--annot', '-a',
                        help='Absolute path for BED file.', default="");

    parser.add_argument('--known_piRNA', 
                        help='Absolute path for piRBase gold piRNA sequences.', required=True);

    args=parser.parse_args();
    
    ## Lets convert bam2bed
    args.path_given = HCGB_files.create_folder(args.path_given)
    
    ## 
    known_piRNA = get_known_piRNA(args.known_piRNA, True)
    
    print(known_piRNA)
    exit()
        
    ## converts BAM file in BED format
    mapping_bed_file = bedtools_caller.convert_bam2bed(args.name, args.input, args.path_given, debug=True)
    
    ## lets split the big file provided
    files_generated = pilfer_module_call(args.path_given, args.name, mapping_bed_file, args.annot, Debug=True)
    
    
    return ()


######
if __name__== "__main__":
    main()