#!/usr/bin/env python3
############################################################
## Jose F. Sanchez, Marta Lopez & Lauro Sumoy             ##
## Copyright (C) 2019-2021 Lauro Sumoy Lab, IGTP, Spain   ##
############################################################
"""
This module performs tRNA analysis with several possible tools: 
MINTMap, ... 

It shows results unifiying formats by using the miRTop nomenclature.
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
from XICRA.modules import database
from XICRA.scripts import generate_DE
from XICRA.scripts import MINTMap_caller

##############################################
def run_tRNA(options):
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
    elif (options.help_tRNA):
        ## information for join reads
        help_XICRA.help_tRNA()
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
    functions.aesthetics_functions.boxymcboxface("tRNA analysis")
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
    print ("+ Software for tRNA analysis selected:")
    print (options.soft_name)
    
    ## get files
    print ('+ Getting files from input folder... ')
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
            pd_samples_retrieved = sampleParser.files.get_files(options, input_dir, "fastq", ("fastq", "fq", "fastq.gz", "fq.gz"), options.debug)
        else:
            print ('+ Mode: join.\n+ Extension: ')
            print ("[_joined.fastq]\n")
            pd_samples_retrieved = sampleParser.files.get_files(options, input_dir, "trim", ['_trim'], options.debug)
    
    ## debug message
    if (Debug):
        print (colored("**DEBUG: pd_samples_retrieve **", 'yellow'))
        print (pd_samples_retrieved)

    ## species
    print ("+ Species provided:", options.species)

    ############################################################
    ## Download tRNA information from MINTmap github site
    ############################################################
    if not (options.database):
        install_path =  os.path.dirname(os.path.realpath(__file__))
        options.database = os.path.join(install_path, "db_files") 
    else:
        options.database = os.path.abspath(options.database)
    
    print ("+ Create folder to store results: ", options.database)
    functions.files_functions.create_folder(options.database)
    
    ## call database module and return tRNA databse generated or updated
    options.tRNA_db = database.tRNA_db(options.database, options.tRNA_db, options.debug)    
    
    ##############################################################
    ## Start the analysis
    ##############################################################
    
    ## generate output folder, if necessary
    if not options.project:
        print ("\n+ Create output folder(s):")
        functions.files_functions.create_folder(outdir)
    
    ## for samples
    outdir_dict = functions.files_functions.outdir_project(outdir, options.project, pd_samples_retrieved, "tRNA", options.debug)
    
    ## optimize threads
    name_list = set(pd_samples_retrieved["new_name"].tolist())
    threads_job = functions.main_functions.optimize_threads(options.threads, len(name_list)) ## threads optimization
    max_workers_int = int(options.threads/threads_job)

    ## to FIX: MINTmap requires to chdir to folder to create results
    max_workers_int = 1

    ## debug message
    if (Debug):
        print (colored("**DEBUG: options.threads " +  str(options.threads) + " **", 'yellow'))
        print (colored("**DEBUG: max_workers " +  str(max_workers_int) + " **", 'yellow'))
        print (colored("**DEBUG: cpu_here " +  str(threads_job) + " **", 'yellow'))

    print ("+ Create a tRNA analysis for each sample retrieved...")    
    
    ## call tRNA_analysis: 
    ## Get user software selection: mintmap, ...
    
    ## dictionary results
    global results_df
    results_df = pd.DataFrame(columns=("name", "soft", "type", "filename"))
    
    # Group dataframe by sample name
    sample_frame = pd_samples_retrieved.groupby(["new_name"])
    
    ## send for each sample
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers_int) as executor:
        commandsSent = { executor.submit(tRNA_analysis, 
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

    print ("\n\n+ tRNA analysis is finished...")
    print ("+ Let's summarize all results...")
    
    ## outdir
    outdir_report = functions.files_functions.create_subfolder("report", outdir)
    expression_folder = functions.files_functions.create_subfolder("tRNA", outdir_report)

    ## debugging messages
    if options.debug:
        print (results_df)
        
    ## merge all parse gtf files created
    print ("+ Summarize tRNA analysis for all samples...")
    
    if 'mintmap' in options.soft_name:
        results_df = results_df.set_index('type')
        
        ## exclusive tRFs
        print ("\n\n+ Parsing exclusive tRNA analysis for all samples...")
        generate_DE.generate_DE(results_df.filter(like="amb", axis=0).set_index('name'), 
                                options.debug, expression_folder,  type_analysis="tRF-amb")
        
        ## amb tRFs
        print ("\n\n+ Parsing ambiguous tRNA analysis for all samples...")
        generate_DE.generate_DE(results_df.filter(like="exc", axis=0).set_index('name'), 
                                options.debug, expression_folder,  type_analysis="tRF-exc")
    else:
        generate_DE.generate_DE(results_df, options.debug, expression_folder,  type_analysis="tRNA")

    print ("\n*************** Finish *******************")
    start_time_partial = functions.time_functions.timestamp(start_time_total)
    print ("\n+ Exiting tRNA module.")
    return()


#########################################
def tRNA_analysis(reads, folder, name, threads, soft_list, species, database, Debug):
    
    ##
    for soft in soft_list:
        if (soft == "mintmap"):
            ## create mintmap
            MINTmap_folder = functions.files_functions.create_subfolder('mintmap', folder)
            code_success = MINTMap_caller.MINTmap_caller(MINTmap_folder, reads, name, threads, species, database, Debug)   
            
            if not code_success:
                print ('** Some error ocurred during MINTmap analysis for sample %s...' %name)
                return ()
            
            ## save results in dataframe
            filename_amb = os.path.join(MINTmap_folder, 'mintmap_parse', name + '_amb.tsv')
            filename_exc = os.path.join(MINTmap_folder, 'mintmap_parse', name + '_exc.tsv')
            
            results_df.loc[len(results_df)] = name, soft, "amb", filename_amb
            results_df.loc[len(results_df)] = name, soft, "exc", filename_exc
    
