#!/usr/bin/env python3
############################################################
## Jose F. Sanchez                                        ##
## Copyright (C) 2019-2020 Lauro Sumoy Lab, IGTP, Spain   ##
############################################################
"""
This module performs miRNA analysis with three possible tools: 
Miraligner, OptimiR and sRNAbench. It shows results unifiying 
formats by using the miRTop nomenclature.
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

from XICRA.scripts import mirtop_caller
from XICRA.scripts import sRNAbench_caller
from XICRA.scripts import optimir_caller
from XICRA.scripts import miraligner_caller


##############################################
def run_miRNA(options):
    """Main function of the module, organizes the miRNA analysis

    First, checks if the files needed to do the mapping (and posterior
    classification of the reads) have been provided by the user:
    - miRNA gff3 annotation, hairpin fasta, mature fasta, miRBase str annotation
    If some is missing it will be downloaded from miRBase. 

    Then, it calls miRNA_analysis() for each sample in parallel.
    Gets user software selection: sRNAbench, optimiR, miraligner.
    Standarize results using miRTop.
    Finally, build final matrix comparing all samples.
    
    :param options: input parameters introduced by the user. See XICRA miRNA -h.

    :returns: None
    """

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
    elif (options.help_miRNA):
        ## information for join reads
        help_XICRA.help_miRNA()
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
    functions.aesthetics_functions.boxymcboxface("miRNA analysis")
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
    print ("+ Software for miRNA analysis selected:")
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
            pd_samples_retrieved = sampleParser.files.get_files(options, input_dir, "fastq", ["fastq", "fq", "fastq.gz", "fq.gz"], options.debug)
        else:
            print ('+ Mode: join.\n+ Extension: ')
            print ("[_joined.fastq]\n")
            pd_samples_retrieved = sampleParser.files.get_files(options, input_dir, "trim", ['_trim'], options.debug)
    
    ## debug message
    if (Debug):
        print (colored("**DEBUG: pd_samples_retrieve **", 'yellow'))
        print (pd_samples_retrieved)

    ## Additional sRNAbench or miRTop options
    
    ## species
    print ("+ Species provided:", options.species)
    
    ############################################################
    ## Download miRNA information: hairpin, mature, str, gff3
    ############################################################
    if not (options.database):
        install_path =  os.path.dirname(os.path.realpath(__file__))
        options.database = os.path.join(install_path, "db_files") 
    else:
        options.database = os.path.abspath(options.database)
    
    print ("+ Create folder to store results: ", options.database)
    functions.files_functions.create_folder(options.database)
    
    ## call database module and return options updated
    options = database.miRNA_db(options)    
    

    
    ##############################################################
    ## Start the analysis
    ##############################################################
    ## generate output folder, if necessary
    if not options.project:
        print ("\n+ Create output folder(s):")
        functions.files_functions.create_folder(outdir)
    
    ## for samples
    outdir_dict = functions.files_functions.outdir_project(outdir, options.project, pd_samples_retrieved, "miRNA", options.debug)
    
    ## optimize threads
    name_list = set(pd_samples_retrieved["new_name"].tolist())
    threads_job = functions.main_functions.optimize_threads(options.threads, len(name_list)) ## threads optimization
    max_workers_int = int(options.threads/threads_job)

    ## debug message
    if (Debug):
        print (colored("**DEBUG: options.threads " +  str(options.threads) + " **", 'yellow'))
        print (colored("**DEBUG: max_workers " +  str(max_workers_int) + " **", 'yellow'))
        print (colored("**DEBUG: cpu_here " +  str(threads_job) + " **", 'yellow'))

    print ("+ Create a miRNA analysis for each sample retrieved...")    
    
    ## call miRNA_analysis: 
    ## Get user software selection: sRNAbench, optimir, ...
    ## Standarize using miRTop
    
    # Group dataframe by sample name
    sample_frame = pd_samples_retrieved.groupby(["new_name"])
    
    ## send for each sample
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers_int) as executor:
        commandsSent = { executor.submit(miRNA_analysis, sorted(cluster["sample"].tolist()), 
                                         outdir_dict[name], name, threads_job, options.miRNA_gff,
                                         options.soft_name, options.matureFasta, options.hairpinFasta, 
                                         options.miRBase_str, options.species, options.database, Debug): name for name, cluster in sample_frame }

        for cmd2 in concurrent.futures.as_completed(commandsSent):
            details = commandsSent[cmd2]
            try:
                data = cmd2.result()
            except Exception as exc:
                print ('***ERROR:')
                print (cmd2)
                print('%r generated an exception: %s' % (details, exc))

    print ("\n\n+ miRNA analysis is finished...")
    print ("+ Let's summarize all results...")
    
    ## outdir
    outdir_report = functions.files_functions.create_subfolder("report", outdir)
    expression_folder = functions.files_functions.create_subfolder("miRNA", outdir_report)

    ## Get results files generated
    
    ## dictionary results
    results_SampleParser = sampleParser.files.get_files(options, input_dir, "miRNA", ["mirtop.tsv"], options.debug)
    results_df = pd.DataFrame(columns=("name", "soft", "filename"))
    results_df['name'] = results_SampleParser['name']
    results_df['soft'] = results_SampleParser['ext']
    results_df['filename'] = results_SampleParser['sample']

    ## debugging messages
    if options.debug:
        print (results_df)
    
    ## merge all parse gtf files created
    print ("+ Summarize miRNA analysis for all samples...")
    generate_DE.generate_DE(results_df, options.debug, expression_folder)

    print ("\n*************** Finish *******************")
    start_time_partial = functions.time_functions.timestamp(start_time_total)
    print ("\n+ Exiting miRNA module.")
    return()

###############
def miRNA_analysis(reads, folder, name, threads, miRNA_gff, soft_list, 
                   matureFasta, hairpinFasta, miRBase_str, species, database, Debug):
    """Passes the sample information to be analyzed by the selected softwares

    Checks the introduced softwares and calls them one by one to analyze the given
    sample. It stores their results in separated folders, which are called as the 
    used software. Finally, it unifies each software output in miRTop format.
    
    :param reads: file with sample reads
    :param folder: output folder
    :param name: sample name
    :param threads: selected threads (by defoult 2)
    :param miRNA_gff: miRNA gff3 annotation file
    :param soft_list: list with the selected softwares 
    :param matureFasta: mature fasta file 
    :param hairpinFasta: hairpin fasta file
    :param miRBase_str: miRBase str annotation 
    :param species: species tag ID. Default: hsa (Homo sapiens)
    :param database: path to store miRNA annotation files downloaded
    :param Debug: display complete log.


    :returns: None
    """
    
    for soft in soft_list:
        if (soft == "sRNAbench"):
            ## create sRNAbench
            sRNAbench_folder = functions.files_functions.create_subfolder('sRNAbench', folder)
            code_success = sRNAbench_caller.sRNAbench_caller(reads, sRNAbench_folder, name, threads, species, Debug) ## Any additional sRNAbench parameter?
                
            if not code_success:
                print ('** miRTop would not be executed for sample %s...' %name)
                return ()
            
            ## create folder for sRNAbench results
            miRTop_folder = functions.files_functions.create_subfolder("sRNAbench_miRTop", folder)
            mirtop_caller.miRTop_caller(sRNAbench_folder, miRTop_folder, name, threads, miRNA_gff, hairpinFasta, 'sRNAbench', species, Debug)
            
        ###
        if (soft == "optimir"):
            ## create OptimiR analysis
            optimir_folder = functions.files_functions.create_subfolder('OptimiR', folder)
            code_success = optimir_caller.optimir_caller(reads, optimir_folder, name, threads, matureFasta, hairpinFasta, miRNA_gff, species, Debug) ## Any additional sRNAbench parameter?
            
            ## create folder for Optimir results
            miRTop_folder = functions.files_functions.create_subfolder("OptimiR_miRTop", folder)
            mirtop_caller.miRTop_caller(optimir_folder, miRTop_folder, name, threads, miRNA_gff, hairpinFasta, 'optimir', species, Debug)
            
        ###
        if (soft == "miraligner"):
            
            ## create OptimiR analysis
            miraligner_folder = functions.files_functions.create_subfolder('miraligner', folder)
            code_success = miraligner_caller.miraligner_caller(reads, miraligner_folder, name, threads, database, species, Debug) 
            
            ## create folder for Optimir results
            miRTop_folder = functions.files_functions.create_subfolder("miraligner_miRTop", folder)
            mirtop_caller.miRTop_caller(miraligner_folder, miRTop_folder, name, threads, miRNA_gff, hairpinFasta, 'seqbuster', species, Debug)
            
