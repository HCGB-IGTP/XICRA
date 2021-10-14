#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez, Marta Lopez & Lauro Sumoy           ##
## Copyright (C) 2019-2021 Lauro Sumoy Lab, IGTP, Spain ##
##########################################################
'''
Calls optimir to create miRNA profile
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
def optimir (reads, outpath, file_name, num_threads, matureFasta, hairpinFasta, miRNA_gff, Debug):
    """Executes the optimiR analyis of the sample.

    Checks if the reads are joined and builds the optimiR command to execute it. 
    
    :param reads: file with sample reads
    :param outpath: output folder
    :param file_name: sample name
    :param num_threads: selected threads (by defoult 2)
    :param matureFasta: mature fasta file 
    :param hairpinFasta: hairpin fasta file
    :param miRNA_gff: miRNA gff3 annotation file
    :param Debug: display complete log.


    :returns: The optimiR output
    """
    optimir_exe = set_config.get_exe("optimir", Debug=Debug)
    sRNAbench_db = os.path.abspath(os.path.join(os.path.dirname(optimir_exe), '..')) ## optimir
    logfile = os.path.join(outpath, 'optimir.log')
    errfile = os.path.join(outpath, 'optimir.err')
    
    if (len(reads) > 1):
        print (colored("** ERROR: Only 1 fastq file is allowed please joined reads before...", 'red'))
        exit()
 
    ## create command  
    cmd = "%s process --fq %s --gff_out -o %s --maturesFasta %s --hairpinsFasta %s --gff3 %s > %s 2> %s" %(
        optimir_exe, reads[0], outpath, matureFasta, hairpinFasta, miRNA_gff, logfile, errfile)
    return(functions.system_call_functions.system_call(cmd))


###############       
def optimir_caller(reads, sample_folder, name, threads, matureFasta, hairpinFasta, miRNA_gff, species, Debug):
    """Passes the sample information to be analyzed by optimiR().

    Checks if the computation has been performed before. If not, it calls
    optimir().
    
    :param reads: file with sample reads
    :param sample_folder: output folder
    :param name: sample name
    :param threads: selected threads (by defoult 2)
    :param matureFasta: mature fasta file 
    :param hairpinFasta: hairpin fasta file
    :param miRNA_gff: miRNA gff3 annotation file
    :param species: species tag ID. Default: hsa (Homo sapiens)
    :param Debug: display complete log.


    :returns: True/False
    """
    # check if previously generated and succeeded
    filename_stamp = sample_folder + '/.success'
    if os.path.isfile(filename_stamp):
        stamp = functions.time_functions.read_time_stamp(filename_stamp)
        print (colored("\tA previous command generated results on: %s [%s -- %s]" %(stamp, name, 'OptimiR'), 'yellow'))
    else:
        # Call sRNAbench
    ## no species option for OptimiR
        code_returned = optimir(reads, sample_folder, name, threads, matureFasta, hairpinFasta, miRNA_gff,  Debug)
        if code_returned:
            functions.time_functions.print_time_stamp(filename_stamp)
        else:
            print ('** Sample %s failed...' %name)
            return(False)
        
    return(True)
