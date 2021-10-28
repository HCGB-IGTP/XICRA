#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez                                      ##
## Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain      ##
##########################################################
'''
Calls FASTQC analysis and parses results generated.
'''
## useful imports
import time
import io
import os
import re
import sys
from sys import argv
from io import open
from termcolor import colored

## import my modules
from HCGB import functions
from XICRA.config import set_config

############
def call_fastqc(path, files, sample, fastqc_bin, threads):    
    ## call system for fastqc sample given
    
    #name = functions.files_functions.create_subfolder(sample, path)
    logFile = path + '/' + sample + '.log'
    
    # create string with files path
    files_string = " ".join(files)
    
    ##print ("+ Calling fastqc for samples...")    
    cmd_fastqc = '%s --extract -t %s -o %s %s > %s 2> %s' %(fastqc_bin, threads, path, files_string, logFile, logFile)
    fastq_code = functions.system_call_functions.system_call( cmd_fastqc )
    
    if not fastq_code:
        print ('** Sample %s failed...' %sample)

    ## send command    
    return (fastq_code)
        
############
def run_module_fastqc(path, files, sample, threads):    
    ## Arguments provided via ARGVs

    ## check if previously done and succeeded
    filename_stamp = path + '/.success'
    if os.path.isfile(filename_stamp):
        stamp = functions.time_functions.read_time_stamp(filename_stamp)
        print (colored("\tA previous command generated results on: %s [%s -- %s]" %(stamp, sample, 'fastqc'), 'yellow'))
    else:
        ## call fastqc
        fastqc_bin = set_config.get_exe('fastqc')
        codeReturn = call_fastqc(path, files, sample, fastqc_bin, threads)

        if codeReturn:
            functions.time_functions.print_time_stamp(filename_stamp)
        
        return ()
