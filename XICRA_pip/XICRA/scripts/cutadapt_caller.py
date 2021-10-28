#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez, Marta Lopez & Lauro Sumoy           ##
## Copyright (C) 2019-2021 Lauro Sumoy Lab, IGTP, Spain ##
##########################################################
'''
Calls cutadapt to trim raw reads
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

#############################################
def caller(list_reads, sample_folder, name, threads, min_read_len, Debug, adapters, extra):
    ## check if previously trimmed and succeeded
    filename_stamp = sample_folder + '/.success'
    if os.path.isfile(filename_stamp):
        stamp = functions.time_functions.read_time_stamp(filename_stamp)
        print (colored("\tA previous command generated results on: %s [%s -- %s]" %(stamp, name, 'cutadapt'), 'yellow'))
    else:
        # Call cutadapt
        cutadapt_exe = set_config.get_exe('cutadapt')
        code_returned = cutadapt(cutadapt_exe, list_reads, sample_folder, name, threads, min_read_len, Debug, adapters, extra)
        if code_returned:
            functions.time_functions.print_time_stamp(filename_stamp)
        else:
            print ('** Sample %s failed...' %name)


#############################################
def cutadapt(cutadapt_exe, reads, path, sample_name, num_threads, min_len_given, Debug, adapters, extra):
    """
    
    :param cutadapt_exe:
    :param reads:
    :param path:
    :param sample_name:
    :param num_threads: 
    :param Debug:
    :param adapters
    :param extra:
    
    :type cutadapt_exe:
    :type reads:
    :type path:
    :type sample_name:
    :type num_threads: 
    :type Debug:
    :type adapters: dictionary
    :type extra: string
    
    """
    logfile = os.path.join(path, sample_name + '.cutadapt.log')
    
    if (len(reads) == 2):
        if not adapters['adapter_a'] or not adapters['adapter_A']:
             print ("** ERROR: Missing adapter information")
             exit()
        
        if extra:
            o_param = os.path.join(path, sample_name + '_temp1_trim_R1.fastq')
            p_param = os.path.join(path, sample_name + '_temp1_trim_R2.fastq')
        else:
            p_param = os.path.join(path, sample_name + '_trim_R2.fastq')
            o_param = os.path.join(path, sample_name + '_trim_R1.fastq')
        
        ## paired-end mode, 15 bps as the min length cutoff
        cmd = '%s -j %s -m %s -a %s -A %s -o %s -p %s %s %s > %s' %(cutadapt_exe,  
                                                                       num_threads, min_len_given, 
                                                                       adapters['adapter_a'], 
                                                                       adapters['adapter_A'], o_param, 
                                                                       p_param, reads[0], reads[1], logfile)
    elif (len(reads) == 1):
        if not adapters['adapter_a']:
             print ("** ERROR: Missing adapter information")
             exit()

        if extra:
            o_param = os.path.join(path, sample_name + '_temp1_trim.fastq')
        else:
            o_param = os.path.join(path, sample_name + '_trim.fastq')
        
        ## single-end mode:
        cmd = '%s -j %s -m %s -a %s -o %s %s > %s' %(cutadapt_exe, num_threads, 
                                                     min_len_given,
                                                     adapters['adapter_a'], 
                                                     o_param, reads[0], logfile)    
    else:
        print ('** Wrong number of files provided for sample: %s...' %sample_name)
        return(False)

    ##
    code = functions.system_call_functions.system_call(cmd)

    ## if additional options, run a second cutadapt command
    ## to ensure this options take effect.
    if (extra):
        if (len(reads) == 2):
            o_param2 = os.path.join(path, sample_name + '_trim_R1.fastq')
            p_param2 = os.path.join(path, sample_name + '_trim_R2.fastq')
        
            ## paired-end mode
            extra_cmd = '%s %s -j %s -m %s -a %s -A %s -o %s -p %s %s %s >> %s' %(cutadapt_exe, 
                                                                                  extra, 
                                                                                  num_threads, 
                                                                                  min_len_given, 
                                                                                  adapters['adapter_a'], 
                                                                                  adapters['adapter_A'],
                                                                                  o_param2, p_param2, 
                                                                                  o_param, p_param, logfile)
        
        elif (len(reads) == 1):
            o_param2 = os.path.join(path, sample_name + '_trim.fastq')
            ## single-end mode:
            extra_cmd = '%s %s -j %s -m %s -a %s -o %s %s >> %s' %(cutadapt_exe, 
                                                                   extra, 
                                                                   num_threads, 
                                                                   min_len_given,
                                                                   adapters['adapter_a'], 
                                                                   o_param2, o_param, logfile)    
        
        code2 = functions.system_call_functions.system_call(extra_cmd)
        
        ## remove: o_param p_param
        if (len(reads) == 2):        
            os.remove(p_param)
       
        os.remove(o_param)
        return (code2)
    
    else:
        return (code)
