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
from HCGB import functions
from XICRA.config import set_config
from XICRA.modules import help_XICRA
from XICRA.scripts import generate_DE, bedtools_caller
from XICRA.scripts import MINTMap_caller
from XICRA.scripts import get_length_distribution

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
    
    functions.aesthetics_functions.pipeline_header('XICRA')
    functions.aesthetics_functions.boxymcboxface("piRNA analysis")
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
    
    ## user software selection
    print ("+ Software for piRNA analysis selected:")
    print (options.soft_name)
    
    ## get files
    print ('+ Getting files from input folder... ')
    print ('+ Check if previously mapping was generated...')
    
    ## retrieve mapping files
    results_SampleParser = sampleParser.files.get_files(options, input_dir, "map", ["Aligned.sortedByCoord.out.bam"], options.debug, bam=True)
    pd_samples_retrieved = results_SampleParser
    del results_SampleParser['dirname']
    del results_SampleParser['ext']
    del results_SampleParser['tag']
    del results_SampleParser['new_name']

    results_SampleParser = results_SampleParser.set_index('name')
    mapping_results = results_SampleParser.to_dict()['sample']

    ## Call mapping if mapping_results is empty
    if not mapping_results:

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
        mapping_outdir_dict = files_functions.outdir_project(outdir, options.project, pd_samples_retrieved, "map", options.debug)
        
        ## debug message
        if (Debug):
            print (colored("**DEBUG: mapping_outdir_dict **", 'yellow'))
            print (mapping_outdir_dict)
    
        # time stamp
        start_time_partial = time_functions.timestamp(start_time_total)
    
        ## optimize threads
        name_list = set(pd_samples_retrieved["new_name"].tolist())
        threads_job = main_functions.optimize_threads(options.threads, len(name_list)) ## threads optimization
        max_workers_int = int(options.threads/threads_job)
        
        ## debug message
        if (Debug):
            print (colored("**DEBUG: options.threads " +  str(options.threads) + " **", 'yellow'))
            print (colored("**DEBUG: max_workers " +  str(max_workers_int) + " **", 'yellow'))
            print (colored("**DEBUG: cpu_here " +  str(threads_job) + " **", 'yellow'))
            
        ##############################################
        ## map Reads
        ##############################################
        (start_time_partial, mapping_results) = map.mapReads_module_STAR(options, pd_samples_retrieved, mapping_outdir_dict, 
                        options.debug, max_workers_int, threads_job, start_time_partial, outdir)
    
        ## debug message
        if (Debug):
             print (colored("**DEBUG: mapping_results **", 'yellow'))
             print (mapping_results)
        
        # time stamp
        start_time_partial = time_functions.timestamp(start_time_partial)

    ############################################################
    ## Download piRNA information: piRBase
    ############################################################
    if not (options.database):
        install_path =  os.path.dirname(os.path.realpath(__file__))
        options.database = os.path.join(install_path, "db_files") 
    else:
        options.database = os.path.abspath(options.database)
    
    print ("+ Create folder to store results: ", options.database)
    functions.files_functions.create_folder(options.database)
    
    ## call database module and return options updated
    options = database.piRNA_db(options)    

    ##############################################################
    ## Start the analysis
    ##############################################################
    ## generate output folder, if necessary
    if not options.project:
        print ("\n+ Create output folder(s):")
        functions.files_functions.create_folder(outdir)
    
    ## for samples
    outdir_dict = functions.files_functions.outdir_project(outdir, options.project, pd_samples_retrieved, "piRNA", options.debug)
    
    ## Mapping is done, process bam files
    print ("+ Create a piRNA analysis for each sample retrieved...")    
    print("+ Intersecting annotation file and mapping file...")
    
    # Group dataframe by sample name
    sample_frame = pd_samples_retrieved.groupby(["new_name"])
    
    ## send for each sample
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers_int) as executor:
        commandsSent = { executor.submit(piRNA_analysis, 
                                         sorted(cluster["sample"].tolist()), 
                                         outdir_dict[name], name, threads_job, 
                                         options.soft_name, options.species, 
                                         options.database, Debug): name for name, cluster in sample_frame }

        for cmd2 in concurrent.futures.as_completed(commandsSent):
            details = commandsSent[cmd2]
            try:
                data = cmd2.result()
            except Exception as exc:
                print ('***ERROR:')
                print (cmd2)
                print('%r generated an exception: %s' % (details, exc))

    print ("\n\n+ piRNA analysis is finished...")
    print ("+ Let's summarize all results...")
    
    ## outdir
    outdir_report = functions.files_functions.create_subfolder("report", outdir)
    expression_folder = functions.files_functions.create_subfolder("piRNA", outdir_report)

    ## merge all parse gtf files created
    print ("+ Summarize piRNA analysis for all samples...")
    
    ## dictionary results
    results_SampleParser = sampleParser.files.get_files(options, input_dir, "piRNA", ["xxx.tsv"], options.debug)
    results_df = pd.DataFrame(columns=("name", "soft", "filename"))
    results_df['name'] = results_SampleParser['name']
    results_df['soft'] = results_SampleParser['ext']
    results_df['filename'] = results_SampleParser['sample']

    ## debugging messages
    if options.debug:
        print (results_df)
    
    print ("\n\n+ Parsing piRNA analysis for all samples...")
    generate_DE.generate_DE(results_df, options.debug, expression_folder,  type_analysis="piRNA")

    print ("\n*************** Finish *******************")
    start_time_partial = functions.time_functions.timestamp(start_time_total)
    
    print ("\n+ Exiting piRNA module.")
    return()


#########################################
def piRNA_analysis(bam_file, folder, name, threads, soft_list, species, database, Debug):
    
    ##
    for soft in soft_list:
        if (soft == "pilfer"):
            ## create pilfer analysis
            pilfer_folder = functions.files_functions.create_subfolder('pilfer', folder)
            code_success = pilfer_caller.pilfer_caller(pilfer_folder, bam_file, name, threads, species, database, Debug)   
            
            if not code_success:
                print ('** Some error ocurred during pilfer analysis for sample %s...' %name)
                return ()
                