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
from XICRA.scripts import RNAbiotype, STAR_caller, multiQC_report, get_length_distribution
from XICRA.other_tools import tools

from HCGB import sampleParser
from HCGB.functions import fasta_functions, time_functions
from HCGB.functions import aesthetics_functions, system_call_functions
from HCGB.functions import files_functions, main_functions

## set to use as a module
## allow multiple software to map
## e.g. bowtie

#########################################
def run_mapping(options):
    print("TODO: Implement this module")

#########################################
def mapReads_module_STAR(options, pd_samples_retrieved, outdir_dict, Debug, 
                    max_workers_int, threads_job, start_time_partial, outdir, multimapping):
    
    """Organizes the mapping of the samples, executed in parallel.

    First, checks if the files needed to do the mapping (and posterior
    classification of the reads) have been provided by the user:
    fasta sequence + annotation or STAR index directory. 

    Then, sends the mapping in parallel for each sample calling mapReads_caller().    
    
    Finally, generate the MultiQC report  of the mapping for each sample.
    
    :param options: input parameters introduced by the user. See XICRA biotype -h.
    :param pd_samples_retrieved: data frame with the information of the samples
    :param outdir_dict: dictionary with the names of the samples and their files
    :param Debug: show extra information of the process
    :param max_workers_int: number of workers for each thread
    :param threads_job: number of threads to do the computation
    :param start_time_partial: time of the beggining of the process
    :param outdir: directory to store the results

    :type Debug: boolean
    :type max_workers_int: int
    :type threads_job: int 
    :type start_time_partial: int
    :type outdir: string


    :returns: None
    """

    # Group dataframe by sample name
    sample_frame = pd_samples_retrieved.groupby(["new_name"])
    
    ## options
    STAR_exe = set_config.get_exe("STAR", Debug=Debug)
    cwd_folder = os.path.abspath("./")
    folder=files_functions.create_subfolder('STAR_files', cwd_folder)

    ## For many samples it will have to load genome index in memory every time.
    ## For a unique sample it will not matter. Take care genome might stay in memory.
    ## Use before loop option LoadAndExit and then:
        ## in loop
        ## Use option LoadAndKeep, set shared memory > 30 Gb
    ## when finished loop Remove memory        
    
    ## check reference
    if (options.fasta):
        print ("+ Genome fasta file provided")
        print ("+ Create genomeDir for later usage...")
        options.fasta = os.path.abspath(options.fasta)
        
        ## create genomeDir
        options.genomeDir = STAR_caller.create_genomeDir(folder, STAR_exe, options.threads, options.fasta, options.limitRAM)
        
    elif (options.genomeDir):
        print ("+ genomeDir provided.")
        options.genomeDir = os.path.abspath(options.genomeDir)
        
    ## remove previous reference genome from memory
    print ("+ Remove genome in memory from previous call... (if any)")
    #STAR_caller.remove_Genome(STAR_exe, options.genomeDir, folder, options.threads)
    
    ## load reference genome
    #STAR_caller.load_Genome(folder, STAR_exe, options.genomeDir, options.threads)

    ## functions.time_functions.timestamp
    start_time_partial = time_functions.timestamp(start_time_partial)
    
    print ("+ Mapping sequencing reads for each sample retrieved...")

    ## send for each sample
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers_int) as executor:
        commandsSent = { executor.submit(mapReads_caller_STAR, sorted(cluster["sample"].tolist()), 
                                         outdir_dict[name], name, threads_job, STAR_exe, 
                                         options.genomeDir, options.limitRAM, Debug, multimapping): name for name, cluster in sample_frame }

        for cmd2 in concurrent.futures.as_completed(commandsSent):
            details = commandsSent[cmd2]
            try:
                data = cmd2.result()
            except Exception as exc:
                print ('***ERROR:')
                print (cmd2)
                print('%r generated an exception: %s' % (details, exc))

    print ("\n\n+ Mapping reads has finished...")
    
    ## functions.time_functions.timestamp
    start_time_partial = time_functions.timestamp(start_time_partial)

    ## remove reference genome from memory
    #STAR_caller.remove_Genome(STAR_exe, options.genomeDir, folder, options.threads)
    
    ## functions.time_functions.timestamp
    start_time_partial = time_functions.timestamp(start_time_partial)

    ## retrieve mapping files
    if options.detached:
        input_dir = outdir
    else:
        input_dir = os.path.abspath(options.input)
    
    results_SampleParser = sampleParser.files.get_files(options, input_dir, "map", ["Aligned.sortedByCoord.out.bam"], options.debug, bam=True)
    del results_SampleParser['dirname']
    del results_SampleParser['ext']
    del results_SampleParser['tag']
    #del results_SampleParser['new_name']

    results_SampleParser = results_SampleParser.set_index('name')
    mapping_results = results_SampleParser.to_dict()['sample']

    ## Create mapping report    
    if (options.skip_report):
        print ("+ No report generation...")
    else:
        print ("\n+ Generating a report using MultiQC module.")
        outdir_report = files_functions.create_subfolder("report", outdir)

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
        
        map_report = files_functions.create_subfolder("STAR", outdir_report)
        multiQC_report.multiQC_module_call(my_outdir_list, "STAR", map_report,"-dd 2")
        print ('\n+ A summary HTML report of each sample is generated in folder: %s' %map_report)

    return(start_time_partial, mapping_results)

#################################
def mapReads_caller_STAR(files, folder, name, threads, STAR_exe, genomeDir, limitRAM_option, Debug, multimapping):
    """Mapping of a given sample with STAR

    First, checks if the trimmed unjoined files exist for the sample and also
    if the calculation has not been done previously. 

    Executes STAR and generate the BAM file of the sample with mapReads script. 
    
    :param files: unjoined trimmed files of the sample
    :param folder: sample folder to store the results
    :param name: sample name
    :param threads: number of threads to do the computation
    :param STAR_exe: check the STAR software is available
    :param genomeDir: path to the genome directory to do the mappig
    :param limitRAM_option: limit RAM bytes to be used in the computation
    :param Debug: show extra information of the process
    :param multimapping: Flag to say whether to use multimapping reads or not
    
    :type folder: string
    :type name: string
    :type threads: int 
    :type start_exe: boolean
    :type genomeDir: string
    :type limitRAM_option: int
    :type Debug: boolean
    :type multimapping: boolean

    :returns: None
    """

    ## check if previously joined and succeeded
    filename_stamp = folder + '/.success'
    if os.path.isfile(filename_stamp):
        stamp = time_functions.read_time_stamp(filename_stamp)
        print (colored("\tA previous command generated results on: %s [%s -- %s]" %(stamp, name, 'STAR'), 'yellow'))
    else:
        ##
        if Debug:
            print ("\n** DEBUG: mapReads_caller options **\n")
            print ("folder: " + folder) 
            print ("name: " + name)
            print ("threads: " + str(threads))
            print ("STAR_exe: " + STAR_exe) 
            print ("genomeDir: " + genomeDir) 
            print ("limitRAM_option: " + str(limitRAM_option))
            print ("files: ")
            print (files)
            
        # Call STAR
        code_returned = STAR_caller.mapReads("LoadAndKeep", files, folder, name, STAR_exe, genomeDir, limitRAM_option, threads, Debug, multimapping)
        
        if (code_returned):
            time_functions.print_time_stamp(filename_stamp)
        else:
            print ("+ Mapping sample %s failed..." %name)
    
    ## return results
    return()
    
