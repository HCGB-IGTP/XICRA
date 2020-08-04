#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez                                        ##
## Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain        ##
##########################################################
"""
Creates Quality check sequence adapters within fastq reads.
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
import cutadapt

## import my modules
from XICRA.scripts import multiQC_report
from XICRA.scripts import fastqc_caller
from XICRA.config import set_config
from XICRA.modules import help_XICRA
from HCGB import sampleParser
from HCGB import functions

##############################################
def run_QC(options):
        ## init time
    start_time_total = time.time()

    ##################################
    ### show help messages if desired    
    ##################################
    if (options.help_format):
        ## help_format option
        sampleParser.help_format()
        exit()
    elif (options.help_project):
        ## information for project
        help_info.project_help()
        exit()
    elif (options.help_multiqc):
        ## information for Multiqc
        multiQC_report.multiqc_help()
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
    
    ## set main header
    functions.aesthetics_functions.pipeline_header()
    functions.aesthetics_functions.boxymcboxface("Quality check")
    print ("--------- Starting Process ---------")
    functions.time_functions.print_time()

    ## absolute path for in & out
    input_dir = os.path.abspath(options.input)
    outdir=""

    ## Project mode as default
    if (options.detached):
        options.project = False
        outdir = os.path.abspath(options.output_folder)
    else:
        options.project = True
        outdir = input_dir        
    
    #fastqc(input_dir, outdir, options, start_time_total)
    functions.aesthetics_functions.boxymcboxface("FASTQC Quality check for samples")
    
    ## get files
    print ('+ Getting files from input folder... ')
    print ('+ Mode: fastq.\n+ Extension: ')
    print ("[ fastq, fq, fastq.gz, fq.gz ]\n")
    pd_samples_retrieved = sampleParser.files.get_files(options, input_dir, "fastq", ("fastq", "fq", "fastq.gz", "fq.gz"), options.debug)

    ## debug message
    if (Debug):
        print (colored("\n**DEBUG: pd_samples_retrieve **", 'yellow'))
        print (pd_samples_retrieved)
        print ("\n")

    ## generate output folder, if necessary
    print ("\n+ Create output folder(s):")
    
    ## if not project, outdir contains the dir to put output
    ## in this case, in some other cases might not occur    
    if not options.project:
        functions.files_functions.create_folder(outdir)
    outdir_dict = functions.files_functions.outdir_project(outdir, options.project, pd_samples_retrieved, "fastqc", options.debug)
    
    print ("+ Checking quality for each sample retrieved...")
    start_time_partial = start_time_total
    
    # Group dataframe by sample name
    sample_frame = pd_samples_retrieved.groupby(["name"])

    ## optimize threads
    name_list = set(pd_samples_retrieved["name"].tolist())
    threads_job = functions.main_functions.optimize_threads(options.threads, len(name_list)) ## threads optimization
    max_workers_int = int(options.threads/threads_job)

    ## debug message
    if (Debug):
        print (colored("**DEBUG: options.threads " +  str(options.threads) + " **", 'yellow'))
        print (colored("**DEBUG: max_workers " +  str(max_workers_int) + " **", 'yellow'))
        print (colored("**DEBUG: cpu_here " +  str(threads_job) + " **", 'yellow'))

    ## send for each sample
    print ("+ Calling fastqc for samples...")    
    with concurrent.futures.ThreadPoolExecutor(max_workers=int(max_workers_int)) as executor:
        commandsSent = { executor.submit(fastqc_caller.run_module_fastqc, outdir_dict[name], sorted( cluster["sample"].tolist() ), name, threads_job): name for name, cluster in sample_frame }
        
        for cmd2 in concurrent.futures.as_completed(commandsSent):
            details = commandsSent[cmd2]
            try:
                data = cmd2.result()
            except Exception as exc:
                print ('***ERROR:')
                print (cmd2)
                print('%r generated an exception: %s' % (details, exc))

    print ("+ FASTQC for samples has finished...")    
    
    ## functions.time_functions.timestamp
    start_time_partial = functions.time_functions.timestamp(start_time_partial)

    if (options.skip_report):
        print ("+ No report generation...")
    else:
        print ("\n+ Generating a report using MultiQC module.")
        outdir_report = functions.files_functions.create_subfolder("report", outdir)

        ## get subdirs generated and call multiQC report module
        givenList = []
        print ("+ Detail information for each sample could be identified in separate folders:")
        
        ## call multiQC report module
        givenList = [ v for v in outdir_dict.values() ]
        my_outdir_list = set(givenList)

        ## debug message
        if (Debug):
            print (colored("\n**DEBUG: my_outdir_list for multiqc report **", 'yellow'))
            print (my_outdir_list)
            print ("\n")
        
        fastqc_report = functions.files_functions.create_subfolder("FASTQC", outdir_report)
        multiQC_report.multiQC_module_call(my_outdir_list, "FASTQC", fastqc_report,"")
        print ('\n+ A summary HTML report of each sample is generated in folder: %s' %fastqc_report)

    print ("\n*************** Finish *******************")
    start_time_partial = functions.time_functions.timestamp(start_time_total)

    print ("+ Exiting qc module.")
    exit()
