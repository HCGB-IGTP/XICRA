#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez, Marta Lopez & Lauro Sumoy           ##
## Copyright (C) 2019-2021 Lauro Sumoy Lab, IGTP, Spain ##
##########################################################
'''
Calls miraligner to create miRNA profile
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
from HCGB.functions import fasta_functions
from HCGB import functions
from XICRA.config import set_config

###############       
def miraligner_caller(reads, sample_folder, name, threads, database, species, Debug):
    """Passes the sample information to be analyzed by miraligner().

    Checks if the computation has been performed before. If not, it calls
    miraligner().
    
    :param reads: file with sample reads
    :param sample_folder: output folder
    :param name: sample name
    :param threads: selected threads (by defoult 2)
    :param database: path to store miRNA annotation files downloaded
    :param species: species tag ID. Default: hsa (Homo sapiens)
    :param Debug: display complete log.


    :returns: True/False
    """
    # check if previously generated and succeeded
    filename_stamp = sample_folder + '/.success'
    if os.path.isfile(filename_stamp):
        stamp = functions.time_functions.read_time_stamp(filename_stamp)
        print (colored("\tA previous command generated results on: %s [%s -- %s]" %(stamp, name, 'miraligner'), 'yellow'))
    else:
        # Call miralinger
        code_returned = miraligner(reads, sample_folder, name, database, species, Debug)
        if code_returned:
            functions.time_functions.print_time_stamp(filename_stamp)
        else:
            print ('** Sample %s failed...' %name)
            return(False)
        
    return(True)

###############       
def miraligner (reads, outpath, file_name, database, species, Debug):
    """Executes the miraligner analyis of the sample

    Checks if the reads are joined and builds the miraligner command to execute it. 
    
    :param reads: file with sample reads
    :param outpath: output folder
    :param file_name: sample name
    :param database: path to store miRNA annotation files downloaded
    :param species: species tag ID. Default: hsa (Homo sapiens)
    :param Debug: display complete log.


    :returns:  The miraligner output
    """    
    miraligner_exe = set_config.get_exe("miraligner", Debug=Debug)
    logfile = os.path.join(outpath, 'miraligner.log')
    
    ## output
    outpath_file = os.path.join(outpath, file_name)
    
    if (len(reads) > 1):
        print (colored("** ERROR: Only 1 fastq file is allowed please joined reads before...", 'red'))
        exit()
    
    ## create tabular information of reads
    tabular_info = os.path.join(outpath, file_name + '-tab.freq.txt')
    fasta_functions.reads2tabular(reads[0], tabular_info)
    
    ## create command 
    java_exe = set_config.get_exe('java', Debug=Debug)
    cmd = '%s -jar %s -db %s -sub 1 -add 3 -trim 3 -s %s -i %s -o %s 2> %s' %(
        java_exe, miraligner_exe, database, species, tabular_info, outpath_file, logfile)
    
    return(functions.system_call_functions.system_call(cmd))
