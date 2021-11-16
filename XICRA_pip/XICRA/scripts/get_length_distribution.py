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
import HCGB.functions.time_functions as HCGB_time
import HCGB.scripts

#######################################################
def get_length_dist(bam_file, GTF_info, name, chr, path_given, debug=False):
    print()
    
    ## convert BAM mapping file into bed
    convert_bam2bed(name, bam_file, path_given, debug)
    
    ## convert GTF to bed
    
    
    ## intersect
    

#######################################################
def convert_GTF2bed(GTF_file, path_given, cpus2use=2, debug=False):
    """
    This functions calls gtf2bed from HCGB to generate a conversion from GTF to BED format.
    """
    
    HCGB_files.create_folder(path_given)
    
    print("+ Split GTF in multiple subsets to speed up process")
            
    ## split GTF file
    split_GTF_files = HCGB.scripts.file_splitter.split_GTF_call(GTF_file, num_files=1, 
                                             name=False, chr_option=True, 
                                             path_given=path_given, debug=debug) 
    # create bed
    print("+ Converting GTF files in BED format:")
    bed_files = {}
    for each_file in split_GTF_files.values():
        ## create name
        each_file_name = HCGB_files.get_file_name(each_file)
        bed_file = os.path.join(path_given, each_file_name + ".bed") ## create a name
        
        print("\t- Converting GTF file: " + each_file)
        print("\t- To BED file: " + each_file)
        print()
        ## convert GTF
        HCGB.scripts.gtf2bed.parse_GTF_call(each_file, bed_file, debug)
        
        ## save
        bed_files[each_file] = bed_file
    
    return (bed_files, split_GTF_files)

#######################################################
def convert_bam2bed(sample, bam_file, path_given, debug):
    """
    This functions calls bedtools to generate a conversion from BAM to BED format.
    """

    bed_file = os.path.join(path_given, HCGB_files.get_file_name(bam_file) + ".bed") ## create a name

    filename_stamp = path_given + '/.convert_bam2bed_success'
    if os.path.isfile(filename_stamp):
        if HCGB_files.is_non_zero_file(bed_file):
            stamp = HCGB_time.read_time_stamp(filename_stamp)
            print (colored("\tA previous command generated results on: %s [%s -- %s]" %(stamp, sample, 'convert_bam2bed'), 'yellow'))
            return (bed_file)
    
    ## execute conversion
    bedtools_exe = set_config.get_exe("bedtools", debug)
    cmd_bedtools = "%s bamtobed -i %s > %s" %(bedtools_exe, bam_file, bed_file)
    
    bed_code = HCGB_sys.system_call(cmd_bedtools, False, True)

    if not bed_code:
        print ("** ERROR: Some error occurred during conversion from BAM to BED... **")
        exit()
    
    ## print time stamp
    HCGB_time.print_time_stamp(filename_stamp)

    return (bed_file)

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
    
    ##
    (bed_files, gtf_files) = convert_GTF2bed(os.path.abspath(args.annot), os.path.abspath(args.path), debug=False)
    
    print("bed_files")
    print(bed_files)
    
    print()
    print("gtf_files")
    print(gtf_files)
    
    exit()
    
    
    ## lets split the big file provided
    files_generated = get_length_dist(os.path.abspath(args.input), os.path.abspath(args.annot), name=args.name, 
              chr=args.annot_splitted, path_given=os.path.abspath(args.path), debug=True)
    return ()

######
if __name__== "__main__":
    main()
    