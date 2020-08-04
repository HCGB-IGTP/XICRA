#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez                                        ##
## Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain        ##
##########################################################
"""
Joins paired-end sequence reads that overlap.
"""
## import useful modules
import os
import sys
import re
import time
from io import open
import shutil
import concurrent.futures
from termcolor import colored

## import my modules
from XICRA.modules import help_XICRA
from XICRA.config import set_config
from HCGB import functions
from HCGB import sampleParser

##############################################
def run_join(options):

    ## init time
    start_time_total = time.time()

    ##################################
    ### show help messages if desired    
    ##################################
    if (options.help_format):
        ## help_format option
        help_XICRA.help_fastq_format()
    elif (options.help_project):
        ## information for project
        help_XICRA.project_help()
        exit()
    elif (options.help_join_reads):
        ## information for join reads
        help_XICRA.help_join_reads()
        exit()
        
    ## debugging messages
    global Debug
    if (options.debug):
        Debug = True
    else:
        Debug = False
        
    ### set as default paired_end mode
    if (options.single_end):
        options.pair = False
    else:
        options.pair = True
    
    functions.aesthetics_functions.pipeline_header()
    functions.aesthetics_functions.boxymcboxface("Join paired-end reads")
    print ("--------- Starting Process ---------")
    functions.time_functions.print_time()

    ## absolute path for in & out
    input_dir = os.path.abspath(options.input)
    outdir=""

    ## set mode: project/detached
    if (options.detached):
        outdir = os.path.abspath(options.output_folder)
        options.project = False
    else:
        options.project = True
        outdir = input_dir        
    
    ## Percentage difference for joining sequences
    if not options.perc_diff:
        options.perc_diff = 0
    
    print ('+ Getting files from input folder... ')
    ## get files
    if options.noTrim:
        print ('+ Mode: fastq.\n+ Extension: ')
        print ("[ fastq, fq, fastq.gz, fq.gz ]\n")
        pd_samples_retrieved = sampleParser.files.get_files(options, input_dir, "fastq", ("fastq", "fq", "fastq.gz", "fq.gz"), options.debug)
    else:
        print ('+ Mode: trim.\n+ Extension: ')
        print ("[ _trim_ ]\n")
        pd_samples_retrieved = sampleParser.files.get_files(options, input_dir, "trim", ['_trim_'], options.debug)
    
    ## debug message
    if (Debug):
        print (colored("**DEBUG: pd_samples_retrieve **", 'yellow'))
        print (pd_samples_retrieved)

    ## generate output folder, if necessary
    print ("\n+ Create output folder(s):")
    if not options.project:
        functions.files_functions.create_folder(outdir)
    ## for samples
    outdir_dict = functions.files_functions.outdir_project(outdir, options.project, pd_samples_retrieved, "join", options.debug)
    
    ## optimize threads
    name_list = set(pd_samples_retrieved["new_name"].tolist())
    threads_job = functions.main_functions.optimize_threads(options.threads, len(name_list)) ## threads optimization
    max_workers_int = int(options.threads/threads_job)

    ## debug message
    if (Debug):
        print (colored("**DEBUG: options.threads " +  str(options.threads) + " **", 'yellow'))
        print (colored("**DEBUG: max_workers " +  str(max_workers_int) + " **", 'yellow'))
        print (colored("**DEBUG: cpu_here " +  str(threads_job) + " **", 'yellow'))

    print ("+ Joining paired-end sequencing reads for each sample retrieved...")    
    
    # Group dataframe by sample name
    sample_frame = pd_samples_retrieved.groupby(["new_name"])
    
    ## send for each sample
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers_int) as executor:
        commandsSent = { executor.submit(fastqjoin_caller, sorted(cluster["sample"].tolist()), 
                                         outdir_dict[name], name, threads_job, options.perc_diff,
                                         Debug): name for name, cluster in sample_frame }

        for cmd2 in concurrent.futures.as_completed(commandsSent):
            details = commandsSent[cmd2]
            try:
                data = cmd2.result()
            except Exception as exc:
                print ('***ERROR:')
                print (cmd2)
                print('%r generated an exception: %s' % (details, exc))

    print ("\n\n+ Joining reads has finished...")
    
    ## TODO: create statistics on joined reads
    ##

    print ("\n*************** Finish *******************")
    start_time_partial = functions.time_functions.timestamp(start_time_total)
    print ("\n+ Exiting join module.")
    return()

#############################################
def fastqjoin_caller(list_reads, sample_folder, name, threads, perc_diff, Debug):
    ## check if previously joined and succeeded
    filename_stamp = sample_folder + '/.success'
    if os.path.isfile(filename_stamp):
        stamp = functions.time_functions.read_time_stamp(filename_stamp)
        print (colored("\tA previous command generated results on: %s [%s -- %s]" %(stamp, name, 'fastqjoin'), 'yellow'))
    else:
        # Call fastqjoin
        fastqjoin_exe = set_config.get_exe('fastqjoin')
        code_returned = fastqjoin(fastqjoin_exe, list_reads, sample_folder, name, threads, perc_diff, Debug)
        if code_returned:
            functions.time_functions.print_time_stamp(filename_stamp)
        else:
            print ('** Sample %s failed...' %name)

#############################################
def fastqjoin (fastqjoin_exe, reads, path, sample_name, num_threads, perc_diff, Debug):
    """
    
    :param fastqjoin_exe:
    :param reads:
    :param path:
    :param sample_name:
    :param num_threads: 
    :param Debug:
    
    :type fastqjoin_exe:
    :type reads:
    :type path:
    :type sample_name:
    :type num_threads: 
    :type Debug:
    
    """
    logfile = os.path.join(path, sample_name + '.fastqjoin.log')
    joined_reads = os.path.join(path, sample_name + '_trim_joined.fastq')
    unjoined_1 = os.path.join(path, sample_name + '_trimmed_unjoin_R1.fastq')
    unjoined_2 = os.path.join(path, sample_name + '_trimmed_unjoin_R2.fastq')
    
    ## check paired-end file
    if (len(reads) == 2):
        cmd = fastqjoin_exe + ' -p %s %s %s -o %s -o %s -o %s > %s' %(perc_diff, reads[0], 
                                                                  reads[1], unjoined_1, unjoined_2, 
                                                                  joined_reads, logfile)
    else:
        print ('** Wrong number of files provided for sample: %s...' %sample_name)
        return(False)

    return(functions.system_call_functions.system_call(cmd))
    