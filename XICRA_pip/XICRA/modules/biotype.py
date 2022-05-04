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
from XICRA.modules import help_XICRA, map
from XICRA.scripts import RNAbiotype, multiQC_report, get_length_distribution
from XICRA.other_tools import tools

from HCGB import sampleParser
from HCGB.functions import fasta_functions, time_functions
from HCGB.functions import aesthetics_functions, system_call_functions
from HCGB.functions import files_functions, main_functions

##############################################
def run_biotype(options):
    """Main function of the module, organizes the biotype analysis

    First, checks if the files needed to do the mapping (and posterior
    classification of the reads) have been provided by the user:
    fasta sequence + annotation or STAR index directory. 

    First, it executes the maping calling map XICRA module.
    Then, it calls RNAbiotype_module_call() to execute featureCounts.
    
    Finally, generate the MultiQC featureCounts report of all samples.
    
    :param options: input parameters introduced by the user. See XICRA biotype -h.

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
    
    aesthetics_functions.pipeline_header('XICRA')
    aesthetics_functions.boxymcboxface("RNA biotype analysis")
    print ("--------- Starting Process ---------")
    time_functions.print_time()

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
    
    ## multimapping:
    if options.no_multiMapping:
        multimapping = False
    else:
        multimapping = True
    
    (start_time_partial, mapping_results) = map.mapReads_module_STAR(options, pd_samples_retrieved, mapping_outdir_dict, 
                    options.debug, max_workers_int, threads_job, start_time_partial, outdir, multimapping)

    ## debug message
    if (Debug):
         print (colored("**DEBUG: mapping_results **", 'yellow'))
         print (mapping_results)
    
    # time stamp
    start_time_partial = time_functions.timestamp(start_time_partial)
    
    ## for samples
    biotype_outdir_dict = files_functions.outdir_project(outdir, options.project, pd_samples_retrieved, "biotype", options.debug)
    
    ## debug message
    if (Debug):
        print (colored("**DEBUG: biotype_outdir_dict **", 'yellow'))
        print (biotype_outdir_dict)
    
    ##############################################
    ## feature counts
    ##############################################

    ## get RNAbiotype information
    RNAbiotype.RNAbiotype_module_call(mapping_results, biotype_outdir_dict, options.annotation, 
                                      options.debug, max_workers_int, threads_job, multimapping, options.stranded)

    # time stamp
    start_time_partial = time_functions.timestamp(start_time_partial)
    
    if (options.skip_report):
        print ("+ No report generation...")
    else:
        print ("\n+ Generating a report using MultiQC module for featureCount analysis.")
        outdir_report = files_functions.create_subfolder("report", outdir)

        ## get subdirs generated and call multiQC report module
        givenList = []
        print ("+ Detail information for each sample could be identified in separate folders:")
        
        ## call multiQC report module
        givenList = [ v for v in biotype_outdir_dict.values() ]
        my_outdir_list = set(givenList)

        ## debug message
        if (Debug):
            print (colored("\n**DEBUG: my_outdir_list for multiqc report **", 'yellow'))
            print (my_outdir_list)
            print ("\n")
        
        featureCount_report = files_functions.create_subfolder("featureCount", outdir_report)
        multiQC_report.multiQC_module_call(my_outdir_list, "featureCount", featureCount_report,"-dd 2")
        print ('\n+ A summary HTML report of each sample is generated in folder: %s' %featureCount_report)

        ### Summarizing RNA biotype information
        biotype_report = files_functions.create_subfolder("biotype", outdir_report)
        single_files_biotype = files_functions.create_subfolder("samples", biotype_report)
        
        ## results
        dict_files = {}
        
        for samples in biotype_outdir_dict:
            featurecount_file = os.path.join(biotype_outdir_dict[samples], 'featureCount.out.tsv')
            if files_functions.is_non_zero_file(featurecount_file):
                dict_files[samples] = featurecount_file
            ## copy pdf
            pdf_plot = main_functions.retrieve_matching_files(biotype_outdir_dict[samples], '.pdf', options.debug)
            if files_functions.is_non_zero_file(pdf_plot[0]):
                shutil.copy(pdf_plot[0], single_files_biotype)
        
        ## collapse all information
        all_data = RNAbiotype.generate_matrix(dict_files)
    
        ## print into excel/csv
        print ('+ Table contains: ', len(all_data), ' entries\n')
        
        ## debugging messages
        if Debug:
            print ("** DEBUG: all_data")
            print (all_data)
        
        ## set abs_csv_outfile to be in report folder
        ## copy or link files for each sample analyzed
        abs_csv_outfile = os.path.join(biotype_report, "summary.csv")
        all_data.to_csv(abs_csv_outfile)
       
        ## create plot: call R and XICRA::stats 
        ## outfile_pdf = os.path.join(biotype_report, "RNAbiotypes_summary.pdf") ## creates automatically the name biotypes-plot.pdf
        
        ## R scripts
        biotype_R_script = tools.R_scripts('plot_RNAbiotype_sum', options.debug)
        rscript = set_config.get_exe("Rscript", options.debug)
        cmd_R_plot = "%s %s -f %s -o %s" %(rscript, biotype_R_script, abs_csv_outfile, biotype_report)
        
        ##
        print ("+ Create summary plot for all samples")
        callCode = system_call_functions.system_call(cmd_R_plot)
    
    ## Create length distribution analysis
    ## get_length_distribution.get_length_dist()
    
    ## Add mirtrace analysis and report
    
    
    
    print ("\n*************** Finish *******************")
    start_time_partial = time_functions.timestamp(start_time_total)
    print ("\n+ Exiting join module.")
    return()

