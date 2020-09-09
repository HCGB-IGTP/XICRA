#!/usr/bin/env python3
############################################################
## Jose F. Sanchez                                        ##
## Copyright (C) 2019-2020 Lauro Sumoy Lab, IGTP, Spain   ##
############################################################
"""
Create RNA biotype analysis using STAR and featureCounts
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
from XICRA.config import set_config
from XICRA.modules import help_XICRA
from XICRA.scripts import RNAbiotype
from XICRA.scripts import mapReads

from HCGB import sampleParser
from HCGB.functions import fasta_functions
from HCGB.functions import time_functions
from HCGB.functions import aesthetics_functions
from HCGB.functions import files_functions

##############################################
def run_biotype(options):

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
    elif (options.help_RNAbiotype):
        ## information for join reads
        RNAbiotype.help_info()
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
    functions.aesthetics_functions.boxymcboxface("RNA biotype analysis")
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
    
    ## get files
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
    mapping_outdir_dict = files_functions.outdir_project(outdir, options.project, pd_samples_retrieved, "map", options.debug)
    
    ## map Reads
    mapReads_module(options, pd_samples_retrieved, mapping_outdir_dict)
    
    ## for samples
    biotype_outdir_dict = files_functions.outdir_project(outdir, options.project, pd_samples_retrieved, "biotype", options.debug)

    ## get RNAbiotype information
    RNAbiotype.RNAbiotype_module_call(samples_dict, biotype_outdir_dict, gtf_file, threads, Debug)

    abs_path_folder = os.path.abspath(argv[1])
    abs_csv_outfile = os.path.abspath(argv[2])
    list_files = os.listdir(abs_path_folder)
    dict_files = {}
    
    for samples in biotype_outdir_dict.items:
        featurecount_file = os.path.join(biotype_outdir_dict[samples], 'featureCount.out.tsv')
        if files_functions.is_non_zero_file(featurecount_file):
            dict_files[l] = featurecount_file
    
    ## collapse all information
    all_data = RNAbiotype.generate_matrix(dict_files)

    ## print into excel/csv
    print ('+ Table contains: ', len(all_data), ' entries\n')
    all_data.to_csv(abs_csv_outfile, quoting=csv.QUOTE_NONNUMERIC)
    
    ## create plot



    
    print ("\n*************** Finish *******************")
    start_time_partial = functions.time_functions.timestamp(start_time_total)
    print ("\n+ Exiting join module.")
    return()


#########################################
def mapReads_module(options, pd_samples_retrieved, outdir_dict):
    
    ## optimize threads
    name_list = set(pd_samples_retrieved["new_name"].tolist())
    threads_job = functions.main_functions.optimize_threads(options.threads, len(name_list)) ## threads optimization
    max_workers_int = int(options.threads/threads_job)

    ## debug message
    if (Debug):
        print (colored("**DEBUG: options.threads " +  str(options.threads) + " **", 'yellow'))
        print (colored("**DEBUG: max_workers " +  str(max_workers_int) + " **", 'yellow'))
        print (colored("**DEBUG: cpu_here " +  str(threads_job) + " **", 'yellow'))

    print ("+ Mapping sequencing reads for each sample retrieved...")    
    
    # Group dataframe by sample name
    sample_frame = pd_samples_retrieved.groupby(["new_name"])
    
        
    ## For many samples it will have to load genome index in memory every time.
    ## For a unique sample it will not matter. Take care genome might stay in memory.
    ## Use before loop option LoadAndExit and then:
        ## in loop
        ## Use option LoadAndKeep, set shared memory > 30 Gb
    ## when finished loop Remove memory        
    ## Send a process for each sample
    
    ## options
    STAR_exe = set_config.get_exe("STAR")
    folder = ""
    
    ## load reference genome
    mapReads.load_Genome(folder, STAR_exe, options.genomeDir, threads_job)
    
    ## send for each sample
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers_int) as executor:
        commandsSent = { executor.submit(mapReads_caller, sorted(cluster["sample"].tolist()), 
                                         outdir_dict[name], name, threads_job, STAR_exe, 
                                         options.genomeDir, options.limitRAM, Debug): name for name, cluster in sample_frame }

        for cmd2 in concurrent.futures.as_completed(commandsSent):
            details = commandsSent[cmd2]
            try:
                data = cmd2.result()
            except Exception as exc:
                print ('***ERROR:')
                print (cmd2)
                print('%r generated an exception: %s' % (details, exc))

    print ("\n\n+ Mapping reads has finished...")
    
    ## remove reference genome from memory
    mapReads.remove_Genome(STAR_exe, options.genomeDir, remove_folder, num_threads)
    
    ## TODO: create statistics on mapped reads

    
    return()

#################################
def mapReads_caller(files, folder, name, threads, STAR_exe, genomeDir, limitRAM_option, Debug):
    ## check if previously joined and succeeded
    filename_stamp = folder + '/.success'
    if os.path.isfile(filename_stamp):
        stamp = time_functions.read_time_stamp(filename_stamp)
        print (colored("\tA previous command generated results on: %s [%s -- %s]" %(stamp, name, 'STAR'), 'yellow'))
    else:
        # Call STAR
        code_returned = mapReads.mapReads("LoadAndKeep", files, folder, name, STAR_exe, genomeDir, limitRAM_option, threads, Debug)
        
        if (code_returned):
             time_functions.print_time_stamp(filename_stamp)
        else:
            print ("+ Mapping sample %s failed..." %name)

    return()
    