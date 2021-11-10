#!/usr/bin/env python3
############################################################
## Jose F. Sanchez                                        ##
## Copyright (C) 2019-2020 Lauro Sumoy Lab, IGTP, Spain   ##
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

from XICRA.config import set_config
import HCGB.functions.aesthetics_functions as HCGB_aes
import HCGB.functions.files_functions as HCGB_files
import HCGB.functions.system_call_functions as HCGB_sys
import HCGB.scripts

#######################################################
def get_length_dist(bam_file, GTF_info, name, chr, path_given, debug=False):
    print()

#######################################################
def convert_GTF2bed(GTF_file, path_given, debug=False):
    print()
    bed_file = os.path.join(path_given, HCGB_files.get_file_name(GTF_file) + ".bed") ## create a name
    data = HCGB.scripts.gtf2bed.parse_GTF(GTF_file, bed_file, debug)
    
    print(data)

#######################################################
def convert_bam2bed(bam_file, path_given, debug):
    """
    This functions calls bedtools to generate a conversion from BAM to BED format.
    """

    bed_file = os.path.join(path_given, HCGB_files.get_file_name(bam_file) + ".bed") ## create a name
    bedtools_exe = set_config.get_exe("bedtools", debug)
    cmd_bedtools = "%s bamtobed -i %s > %s" %(bedtools_exe, bam_file, bed_file)
    
    bed_code = HCGB_sys.system_call(cmd_bedtools, False, True)

    if not bed_code:
        print ("** ERROR: Some error occurred during genomeDir creation... **")
        exit()
    
    return (bed_code)

#######################################################
def main():
    ## this code runs when call as a single script
    parser=argparse.ArgumentParser(description='''Get length distribution for biotypes mapped in BAM file
    Note:
    It takes into account no genes or transcript are broken
    It is also possible to split according to chromosome (one gtf/chromosome)
    
    ''');
    
    parser.add_argument('--input', '-i', help='Input BAM file', required=True);

    parser.add_argument('--name',
                        help='Name of the sample. Default: use filename provided.', default="");
    
    parser.add_argument('--path',
                        help='Path to save results generated. Default: use path from filename provided.', default="");

    parser.add_argument('--annot',
                        help='Absolute path for GTF file.', default="");

    parser.add_argument('--annot_splitted', action="store_true",
                        help='--annot is a dir containing splitted GTF files for each chromosome.');
    
    args=parser.parse_args();
    
    ## lets split the big file provided
    files_generated = get_length_dist(os.path.abspath(args.input), os.path.abspath(args.annot), name=args.name, 
              chr=args.annot_splitted, path_given=args.path, debug=True)
    
    return ()

######
if __name__== "__main__":
    main()
    