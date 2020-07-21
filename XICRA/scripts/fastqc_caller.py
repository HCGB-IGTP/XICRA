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

## import my modules
from XICRA.scripts import functions
from XICRA.config import set_config

############
def call_fastqc(path, file1, file2, sample, fastqc_bin, threads):    
    ## call system for fastqc sample given
    #name = functions.create_subfolder(sample, path)
    logFile = path + '/' + sample + '.log'
    
    if os.path.isfile(logFile):
        return ('OK')
    
    ##print ("+ Calling fastqc for samples...")    
    cmd_fastqc = '%s --extract -t %s -o %s %s %s > %s 2> %s' %(fastqc_bin, threads, path, file1, file2, logFile, logFile)
    ## send command    
    return (functions.system_call( cmd_fastqc ))
        
############
def run_module_fastqc(path, files, sample, threads):    
    ## Arguments provided via ARGVs
    fastqc_bin = set_config.get_exe('fastqc')

    ## check if paired end
    if (len(files) < 2):
        print ('+ No implementation yet for single-end. Sorry.')
        exit()

    codeReturn = call_fastqc(path, files[0], files[1], sample, fastqc_bin, threads)
    if codeReturn == 'FAIL':
        exit()
    path_to_sample = path + '/' + sample
    return path_to_sample
