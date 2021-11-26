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
def bam2sam(sample, input_file, path_given, ncpu, options, sam2bam=False, Debug=False):
    """
    This functions calls samtools to generate a conversion from BAM to SAM format.
    """
    ## create sam_file
    outfile_tmp = HCGB_files.get_path_name(input_file, path_given, sample, debug=Debug)
    
    if sam2bam:
        outfile = outfile_tmp.split('.sam')[0] + '.bam'
    else:
        outfile = outfile_tmp.split('.bam')[0] + '.sam'
    
    if not path_given:
        path_given = os.path.dirname(bam_file)
    
    filename_stamp = path_given + '/.convert_bam2sam_success'
    if os.path.isfile(filename_stamp):
        if HCGB_files.is_non_zero_file(outfile):
            stamp = HCGB_time.read_time_stamp(filename_stamp)
            print (colored("\tA previous command generated results on: %s [%s -- %s]" %(stamp, sample, 'samtools view conversion'), 'yellow'))
            return (outfile)

    ## Create call 
    samtools_exe = set_config.get_exe("samtools", Debug)
    cmd_samtools = "%s view --threads %s -o %s" %(samtools_exe, ncpu, outfile)
    
    if (options):
        cmd_samtools = cmd_samtools + " " + str(options)
    
    if (sam2bam):
        cmd_samtools = cmd_samtools + " -b"
    
    cmd_samtools = cmd_samtools + " " + input_file
    samtools_code = HCGB_sys.system_call(cmd_samtools, False, True)

    if samtools_code:
        ## print time stamp
        HCGB_time.print_time_stamp(filename_stamp)
    
        return (outfile)

