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
from termcolor import colored

from XICRA.config import set_config
import HCGB.functions.aesthetics_functions as HCGB_aes
import HCGB.functions.files_functions as HCGB_files
import HCGB.functions.system_call_functions as HCGB_sys
import HCGB.functions.time_functions as HCGB_time
import HCGB.scripts

import pybedtools

#######################################################
def get_length_dist(bam_file, GTF_info, name, chr, path_given, folder_GTF, debug=False):
    
    ## debugging messages
    if debug:
        HCGB_aes.debug_message("******************")
        HCGB_aes.debug_message("get_length_dist function call: ", color="yellow")
        print()
        HCGB_aes.debug_message("bam_file: " + bam_file, color="yellow")
        HCGB_aes.debug_message("GTF_info: " + GTF_info, color="yellow")
        HCGB_aes.debug_message("name: " + name, color="yellow")
        HCGB_aes.debug_message("path_given: " + path_given, color="yellow")
        HCGB_aes.debug_message("folder_GTF: " + folder_GTF, color="yellow")
        HCGB_aes.debug_message("chr: " + str(chr), color="yellow")
        
    
    ## converts BAM file in BED format
    sample_bed_file = convert_bam2bed(name, bam_file, path_given, debug=debug)
    
    ## convert annotation in GTF to BED
    if chr: 
        ## GTF_info is already splitted:
        (bed_files, gtf_files) = convert_GTF2bed(GTF_info, folder_GTF, num_files=0, debug=debug)
    else:    
        (bed_files, gtf_files) = convert_GTF2bed(GTF_info, folder_GTF, debug=debug)

    ## debugging messages
    if debug:
        HCGB_aes.debug_message("bed_files", color="yellow")
        print(bed_files)        
        print()
        HCGB_aes.debug_message("gtf_files", color="yellow")
        print(gtf_files)
        
    ## Loop for each reference sequence BED file
    ## and intersect with mapping BAM->BED file
    for RefSeq_bed_file in bed_files.values():
        print(RefSeq_bed_file)

        ## intersect
    
    exit()

    
#######################################################
def convert_GTF2bed(GTF_file, path_given, cpus2use=2, num_files=1, debug=False):
    """
    This functions calls gtf2bed from HCGB to generate a conversion from GTF to BED format.
    """
    
    HCGB_files.create_folder(path_given)
    
    split_GTF_files = {}
    if num_files==0:
        ## GTF file is already splitted
        split_GTF_files['GTF_file'] = os.path.abspath(GTF_file)
    else:
        
        print("+ Split GTF in multiple subsets to speed up process")
                
        ## split GTF file
        split_GTF_files = HCGB.scripts.file_splitter.split_file_call(GTF_file, num_files=1, 
                                                 name=False, chr_option=True, 
                                                 in_format="GTF",
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
    
    It converts, sorts and collapse information retaining counts for each feature in bed format
    """

    bed_file = os.path.join(path_given, HCGB_files.get_file_name(bam_file) + ".bed") ## create a name
    bed_file_tmp = bed_file + '_tmp' 
    
    filename_stamp = path_given + '/.convert_bam2bed_success'
    if os.path.isfile(filename_stamp):
        if HCGB_files.is_non_zero_file(bed_file):
            stamp = HCGB_time.read_time_stamp(filename_stamp)
            print (colored("\tA previous command generated results on: %s [%s -- %s]" %(stamp, sample, 'convert_bam2bed'), 'yellow'))
            return (bed_file)
    
    ## execute conversion and count reads mapping in exact coordinates
    ## bedtools bamtobed -i bam_file > bed_file 
    ## bedtools sort -chrThenSizeA -i bed_file > bed_sort.file
    ## bedtools groupby -i bed_sort.file -o count -g 1,2,3 -c 4 > counts.bed 
    
    bedtools_exe = set_config.get_exe("bedtools", debug)
    cmd_bedtools = "%s bamtobed -i %s | %s groupby -o count -g 1,2,3 -c 4 > %s" %(bedtools_exe, bam_file, bedtools_exe, bed_file_tmp)
    bed_code = HCGB_sys.system_call(cmd_bedtools, False, True)

    if bed_code:
        cmd_bedtools2 = "%s sort -chrThenSizeA -i %s | %s groupby -o count -g 1,2,3 -c 4 > %s" %(bedtools_exe, bed_file_tmp, 
                                                                                                                bedtools_exe, bed_file)
        bed_code2 = HCGB_sys.system_call(cmd_bedtools2, False, True)
    
    ## -----------------------------------------------
    ## Pybedtools
    ## -----------------------------------------------
    ## It might be possible to use pybedtools. We need to load bam, it is not possible to use string to absolute path
    #bed_info = pybedtools.bedtool.BedTool.bam_to_bed(bam_file)
    
    if not bed_code or not bed_code2:
        print ("** ERROR: Some error occurred during conversion from BAM to BED... **")
        exit()
    
    
    ## print time stamp
    HCGB_time.print_time_stamp(filename_stamp)

    ## remove tmp files
    os.remove(bed_file_tmp)

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

    parser.add_argument('--name', '-n',
                        help='Name of the sample. Default: use filename provided.', default="");
    
    parser.add_argument('--path', '-p',
                        help='Path to save results generated. Default: use path from filename provided.', default="");

    parser.add_argument('--path_GTF',
                        help='Path to store results generated from split GTF file. Default: use path provided previously.', default="");

    parser.add_argument('--annot', '-a',
                        help='Absolute path for GTF file.', default="");

    parser.add_argument('--annot_splitted', '-as', action="store_true",
                        help='--annot is a dir containing splitted GTF files for each chromosome.');
    
    args=parser.parse_args();
    
    ## lets split the big file provided
    files_generated = get_length_dist(os.path.abspath(args.input), os.path.abspath(args.annot), name=args.name, 
              chr=args.annot_splitted, path_given=os.path.abspath(args.path), folder_GTF=os.path.abspath(args.path_GTF), 
              debug=False)
    return ()

######
if __name__== "__main__":
    main()
    