#!/usr/bin/env python3
############################################################
## Jose F. Sanchez, Marta Lopez & Lauro Sumoy             ##
## Copyright (C) 2019-2021 Lauro Sumoy Lab, IGTP, Spain   ##
############################################################
"""
This module performs piRNA analysis with several possible tools: 

Show a description of piRNA molecules here: http://pirnabank.ibab.ac.in/about.html
"""
## import useful modules
import os
import sys
import re
import time
from io import open
import shutil
import concurrent.futures
import pandas as pd
from termcolor import colored

## import my modules
from HCGB import sampleParser

import HCGB.functions.files_functions as HCGB_files
import HCGB.functions.aesthetics_functions as HCGB_aes
import HCGB.functions.time_functions as HCGB_time
import HCGB.functions.main_functions as HCGB_main

from XICRA.config import set_config
from XICRA.modules import help_XICRA, map, database
from XICRA.scripts import pilfer_caller

##############################################
def run_piRNA(options):
    '''
    '''
    
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
    
    elif (options.help_piRNA):
        ## information for join reads
        help_XICRA.help_piRNA()
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
    
    HCGB_aes.pipeline_header('XICRA')
    HCGB_aes.boxymcboxface("piRNA analysis")
    print ("--------- Starting Process ---------")
    HCGB_time.print_time()

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
    
    ## user software selection
    print ("+ Software for piRNA analysis selected:")
    print (options.soft_name)
    
    ## get files
    print ('+ Getting files from input folder... ')
    print ('+ Check if previously mapping was generated...')
    
    ## retrieve mapping files
    results_SampleParser = sampleParser.files.get_files(options, input_dir, "map", ["bam"], options.debug, bam=True)
    pd_samples_retrieved = results_SampleParser
    del results_SampleParser['dirname']
    del results_SampleParser['ext']
    del results_SampleParser['tag']

    results_SampleParser = results_SampleParser.set_index('name')
    mapping_results = results_SampleParser.to_dict()['sample']
    
    ##
    if Debug:
        HCGB_aes.debug_message("results_SampleParser", "yellow")
        print(results_SampleParser)
        
        HCGB_aes.debug_message("mapping_results", "yellow")
        print(mapping_results)        

    ## Call mapping if mapping_results is empty
    if not mapping_results:
        
        if Debug:
            HCGB_aes.debug_message("********************")
            HCGB_aes.debug_message("mapping samples")

        ## get files: use joined reads if paired-end data
        if options.pair:
            options.pair = False ## set paired-end to false for further prepocessing
            if options.noTrim:
                print ('+ Mode: fastq.\n+ Extension: ')
                print ("[ fastq, fq, fastq.gz, fq.gz ]\n")
                pd_samples_retrieved = sampleParser.files.get_files(options, input_dir, "fastq", ["fastq", "fq", "fastq.gz", "fq.gz"], options.debug)
            else:
                print ('+ Mode: join.\n+ Extension: ')
                print ("[_joined.fastq]\n")
                pd_samples_retrieved = sampleParser.files.get_files(options, input_dir, "join", ['_joined.fastq'], options.debug)
        else:
            if options.noTrim:
                print ('+ Mode: fastq.\n+ Extension: ')
                print ("[ fastq, fq, fastq.gz, fq.gz ]\n")
                pd_samples_retrieved = sampleParser.files.get_files(options, input_dir, "fastq", ["fastq", "fq", "fastq.gz", "fq.gz"], options.debug)
            else:
                print ('+ Mode: join.\n+ Extension: ')
                print ("[_joined.fastq]\n")
                pd_samples_retrieved = sampleParser.files.get_files(options, input_dir, "trim", ['_trim'], options.debug)
        
        ## debug message
        if (Debug):
            print (colored("**DEBUG: pd_samples_retrieve **", 'yellow'))
            print (pd_samples_retrieved)
    
        ## generate output folder, if necessary
        print ("\n+ Create output folder(s):")
        if not options.project:
            files_functions.create_folder(outdir)
    
        ## for samples
        mapping_outdir_dict = HCGB_files.outdir_project(outdir, options.project, pd_samples_retrieved, "map", options.debug)
        
        ## debug message
        if (Debug):
            print (colored("**DEBUG: mapping_outdir_dict **", 'yellow'))
            print (mapping_outdir_dict)
    
        # time stamp
        start_time_partial = HCGB_time.timestamp(start_time_total)
    
        ## optimize threads
        name_list = set(pd_samples_retrieved["new_name"].tolist())
        threads_job = HCGB_main.optimize_threads(options.threads, len(name_list)) ## threads optimization
        max_workers_int = int(options.threads/threads_job)
        
        ## debug message
        if (Debug):
            print (colored("**DEBUG: options.threads " +  str(options.threads) + " **", 'yellow'))
            print (colored("**DEBUG: max_workers " +  str(max_workers_int) + " **", 'yellow'))
            print (colored("**DEBUG: cpu_here " +  str(threads_job) + " **", 'yellow'))
            
        ##############################################
        ## map Reads
        ##############################################
        
        ## TODO: Fix options require options.fasta, provided it or options.genomeDir
        (start_time_partial, mapping_results) = map.mapReads_module_STAR(options, pd_samples_retrieved, mapping_outdir_dict, 
                        options.debug, max_workers_int, threads_job, start_time_partial, outdir)
    
        ## debug message
        if (Debug):
             print (colored("**DEBUG: mapping_results **", 'yellow'))
             print (mapping_results)
        
        # time stamp
        start_time_partial = HCGB_time.timestamp(start_time_partial)

    ############################################################
    ## Download piRNA information: piRBase
    ############################################################
    if not (options.database):
        install_path =  os.path.dirname(os.path.realpath(__file__))
        options.database = os.path.join(install_path, "db_files") 
    else:
        options.database = os.path.abspath(options.database)
    
    print ("+ Create folder to store results: ", options.database)
    HCGB_files.create_folder(options.database)
    
    ## call database module and return information updated
    db_info = database.piRNA_info(options.database, options.species, Debug)

    ##############################################################
    ## Start the analysis
    ##############################################################
    ## generate output folder, if necessary
    if not options.project:
        print ("\n+ Create output folder(s):")
        HCGB_files.create_folder(outdir)
    
    ## for samples
    outdir_dict = HCGB_files.outdir_project(outdir, options.project, pd_samples_retrieved, "piRNA", options.debug)
    
    ## Mapping is done, process bam files
    print ("+ Create a piRNA analysis for each sample retrieved...")    
    print("+ Intersecting annotation file and mapping file...")
    
    ## optimize threads
    name_list = set(pd_samples_retrieved["new_name"].tolist())
    #threads_job = HCGB_main.optimize_threads(options.threads, len(name_list)) ## threads optimization
    threads_job = options.threads ## So far all CPUs provided for each sample at a time
    max_workers_int = 1

    ## debug message
    if (Debug):
        print (colored("**DEBUG: options.threads " +  str(options.threads) + " **", 'yellow'))
        print (colored("**DEBUG: max_workers " +  str(max_workers_int) + " **", 'yellow'))
        print (colored("**DEBUG: cpu_here " +  str(threads_job) + " **", 'yellow'))
    # Group dataframe by sample name
    sample_frame = pd_samples_retrieved.groupby(["new_name"])
    
    ## for each sample do piRNA analysis
    for name, cluster in sample_frame:
        if Debug:
            print(name)
            print(cluster)
            
        ## send for each sample        
        piRNA_analysis(sorted(cluster["sample"].tolist())[0], outdir_dict[name], name, threads_job, 
                                         options.soft_name, options.species, options.database, Debug)

    ## outdir
    outdir_report = HCGB_files.create_subfolder("report", outdir)
    expression_folder = HCGB_files.create_subfolder("piRNA", outdir_report)
    
    ## TODO: control if any other software implemented
    ## if pilfer
    pilfer_folder = HCGB_files.create_subfolder("pilfer", expression_folder)
    
    ## get results
    results_pilfer = sampleParser.files.get_files(options, input_dir, "piRNA", ["pilfer_clusters.bed"], options.debug, bam=False)
    del results_pilfer['dirname']
    del results_pilfer['ext']
    del results_pilfer['tag']
    results_pilfer = results_pilfer.set_index('name')
    pilfer_results = results_pilfer.to_dict()['sample']

    ## merge clusters generated
    print ("+ Summarize piRNA analysis for all samples...")
    pilfer_caller.pilfer_merge_samples_call(pilfer_results, pilfer_folder, options.debug)
        
    print ("\n\n+ piRNA analysis is finished...")
    print ("+ Let's summarize all results...")
    
    exit()
    
    print ("\n*************** Finish *******************")
    start_time_partial = HCGB_time.timestamp(start_time_total)
    
    print ("\n+ Exiting piRNA module.")
    return()

#########################################
def piRNA_analysis(bam_file, folder, name, threads, soft_list, species, database, Debug):
    
    if Debug:
        HCGB_aes.debug_message("BAM_FILE: " + bam_file, "yellow")
        HCGB_aes.debug_message("folder: " + folder, "yellow")
        HCGB_aes.debug_message("name: " + name, "yellow")
        HCGB_aes.debug_message("threads: " + str(threads), "yellow")
        HCGB_aes.debug_message("soft_list: ", "yellow")
        print(soft_list)
        HCGB_aes.debug_message("species: " + species, "yellow")
    
    ##
    for soft in soft_list:
        if (soft == "pilfer"):
            ## create pilfer analysis
            pilfer_folder = HCGB_files.create_subfolder('pilfer', folder)
            code_success = pilfer_caller.pilfer_module_call(pilfer_folder, name, bam_file, database, threads, species, Debug)   
            
            if not code_success:
                print ('** Some error ocurred during pilfer analysis for sample %s...' %name)
                return ()
                