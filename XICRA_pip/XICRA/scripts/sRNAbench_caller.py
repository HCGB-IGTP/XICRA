#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez, Marta Lopez & Lauro Sumoy           ##
## Copyright (C) 2019-2021 Lauro Sumoy Lab, IGTP, Spain ##
##########################################################
'''
Calls sRNAbench to create miRNA profile
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


###############
def sRNAbench_caller(reads, sample_folder, name, threads, species, Debug):
    """Passes the sample information to be analyzed by sRNAbench()

    Checks if the computation has been performed before. If not, it calls
    sRNAbench().
    
    :param reads: file with sample reads
    :param sample_folder: output folder
    :param name: sample name
    :param threads: selected threads (by defoult 2)
    :param species: species tag ID. Default: hsa (Homo sapiens)
    :param Debug: display complete log.


    :returns: True/False
    """
    # check if previously generated and succeeded
    filename_stamp = sample_folder + '/.success'
    if os.path.isfile(filename_stamp):
        stamp = functions.time_functions.read_time_stamp(filename_stamp)
        print (colored("\tA previous command generated results on: %s [%s -- %s]" %(stamp, name, 'sRNAbench'), 'yellow'))
    else:
        # Call sRNAbench
        code_returned = sRNAbench(reads, sample_folder, name, threads, species, Debug)
        if code_returned:
            functions.time_functions.print_time_stamp(filename_stamp)
        else:
            print ('** Sample %s failed...' %name)
            return(False)
        
    return(True)

###############       
def sRNAbench (reads, outpath, file_name, num_threads, species, Debug):
    """Executes the sRNAbench analyis of the sample.

    Checks if the reads are joined and builds the sRNAbench command to execute it. 
    
    :param reads: file with sample reads
    :param outpath: output folder
    :param file_name: sample name
    :param num_threads: selected threads (by defoult 2)
    :param species: species tag ID. Default: hsa (Homo sapiens)
    :param Debug: display complete log.


    :returns: The sRNAbench output
    """
    sRNAbench_exe = set_config.get_exe("sRNAbench", Debug=Debug)
    
    ## set as option
    ## sRNAbench.jar in exec/ folder within sRNAtoolboxDB
    ## sRNAbench_db = os.path.abspath(os.path.join(os.path.dirname(sRNAbench_exe), '..')) ## sRNAtoolboxDB
    
    ## sRNAbench.jar linked in bin/. sRNAtoolboxDB in same folder 
    sRNAbench_db = os.path.abspath(os.path.join(os.path.dirname(sRNAbench_exe), 'sRNAtoolboxDB')) ## sRNAtoolboxDB
    
    logfile = os.path.join(outpath, 'sRNAbench.log')
    
    if (len(reads) > 1):
        print (colored("** ERROR: Only 1 fastq file is allowed please joined reads before...", 'red'))
        exit()
    
    ## create command    
    java_exe = set_config.get_exe('java', Debug=Debug)
    cmd = '%s -jar %s dbPath=%s input=%s output=%s' %(java_exe, sRNAbench_exe, sRNAbench_db, reads[0], outpath)  
    cmd = cmd + ' microRNA=%s isoMiR=true plotLibs=true graphics=true' %species 
    cmd = cmd + ' plotMiR=true bedGraphMode=true writeGenomeDist=true'
    cmd = cmd + ' chromosomeLevel=true chrMappingByLength=true > ' + logfile 
    
    return(functions.system_call_functions.system_call(cmd))
