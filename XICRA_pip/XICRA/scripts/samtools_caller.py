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

###################################################
def bam2sam(sample, bam_file, path_given, ncpu, header=True, Debug=False):
    """
    This functions calls samtools to generate a conversion from BAM to SAM format.
    """
    ## create sam_file
    sam_file_tmp = HCGB_files.get_path_name(bam_file, path_given, os.path.basename(bam_file), debug=Debug)
    sam_file = sam_file_tmp.split('.bam')[0] + '.sam'
    
    if not path_given:
        path_given = os.path.dirname(bam_file)
    
    filename_stamp = path_given + '/.convert_bam2sam_success'
    if os.path.isfile(filename_stamp):
        if HCGB_files.is_non_zero_file(sam_file):
            stamp = HCGB_time.read_time_stamp(filename_stamp)
            if HCGB_files.is_non_zero_file(sam_file):
                print (colored("\tA previous command generated results on: %s [%s -- %s]" %(stamp, sample, 'convert_bam2sam'), 'yellow'))
                return (sam_file)

    ## Create call 
    samtools_exe = set_config.get_exe("samtools", Debug)
    if (header):
        cmd_samtools = "%s view --threads %s -h %s > %s" %(samtools_exe, ncpu, bam_file, sam_file) ## include header
    else:
        cmd_samtools = "%s view --threads %s %s > %s" %(samtools_exe, ncpu, bam_file, sam_file)
        
    samtools_code = HCGB_sys.system_call(cmd_samtools, False, True)

    if samtools_code:
        ## print time stamp
        HCGB_time.print_time_stamp(filename_stamp)
    
        return (sam_file)
