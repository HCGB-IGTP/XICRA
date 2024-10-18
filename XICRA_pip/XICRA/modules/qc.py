#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez, Marta Lopez and Lauro Sumoy         ##
## Copyright (C) 2019-2024 Lauro Sumoy Lab, IGTP, Spain ##
##########################################################
"""
Creates Quality check sequence adapters within fastq reads.
"""
## import useful modules
import os
import time
import concurrent.futures
from termcolor import colored


## import my modules
from XICRA import __version__ as pipeline_version
from XICRA.scripts import multiQC_report
from XICRA.scripts import fastqc_caller
from XICRA.modules import help_XICRA
from HCGB import sampleParser
import HCGB.functions.aesthetics_functions as HCGB_aes
import HCGB.functions.time_functions as HCGB_time
import HCGB.functions.files_functions as HCGB_files
import HCGB.functions.info_functions as HCGB_info
import HCGB.functions.main_functions as HCGB_main

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
        help_XICRA.project_help()
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
    HCGB_aes.pipeline_header('XICRA')
    HCGB_aes.boxymcboxface("Quality check")
    print ("--------- Starting Process ---------")
    HCGB_time.print_time()

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
    
    HCGB_aes.boxymcboxface("FASTQC Quality check for samples")
    
    ## get files
    print ('+ Getting files from input folder... ')
    print ('+ Mode: fastq.\n+ Extension: ')
    print ("[ fastq, fq, fastq.gz, fq.gz ]\n")
    pd_samples_retrieved = sampleParser.files.get_files(options, input_dir, "fastq", ["fastq", "fq", "fastq.gz", "fq.gz"], options.debug)

    ## create FASTQC call
    outdir_dict = fastqc(pd_samples_retrieved, outdir, options, "", start_time_total, Debug)

    print ("\n*************** Finish *******************")
    HCGB_time.timestamp(start_time_total)

    ## samples information dictionary
    samples_info = {}
    samples_frame = pd_samples_retrieved.groupby('new_name')
    for name_tuple, grouped in samples_frame:
        #name = name_tuple[0]
        samples_info[name_tuple] = grouped['sample'].to_list()
    
    ## dump information and parameters
    info_dir = HCGB_files.create_subfolder("info", outdir)
    print("+ Dumping information and parameters")
    runInfo = { "module":"qc", "time":time.time(),
                "XICRA version":pipeline_version,
                'sample_info': samples_info,
                'outdir_dict': outdir_dict}
    
    HCGB_info.dump_info_run(info_dir, "qc", options, runInfo, options.debug)
    
    ## dump conda details
    HCGB_info.dump_info_conda(info_dir, "qc", package_name="XICRA", debug=options.debug)


    print ("+ Exiting qc module.")
    exit()

#######################
def fastqc(pd_samples_retrieved, outdir, options, name_analysis, time_stamp, Debug):
    '''
    This is a main function to prepare data to call FASTQC.
    
    :param pd_samples_retrieved
    :param outdir
    :param options
    :param name_analysis
    :param Debug
    
    :type pd_samples_retrieved
    :type outdir
    :type options
    :type name_analysis
    :type Debug
    
    '''
    
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
        HCGB_files.create_folder(outdir)
        
    ## folder name
    if (name_analysis):
        fold_name = "fastqc_" + name_analysis
    else:
        fold_name = "fastqc"

    ## create output dirs for each sample    
    outdir_dict = HCGB_files.outdir_project(outdir, options.project, pd_samples_retrieved, fold_name, options.debug, groupby_col="new_name")
    
    print ("+ Checking quality for each sample retrieved...")
    start_time_partial = time_stamp
    
    # Group dataframe by sample name
    sample_frame = pd_samples_retrieved.groupby(["name"])

    ## optimize threads
    name_list = set(pd_samples_retrieved["name"].tolist())
    threads_job = HCGB_main.optimize_threads(options.threads, len(name_list)) ## threads optimization
    max_workers_int = int(options.threads/threads_job)

    ## debug message
    if (Debug):
        HCGB_aes.debug_message("options.threads: " + str(options.threads), "yellow")
        HCGB_aes.debug_message("max_workers: " + str(max_workers_int), "yellow")
        HCGB_aes.debug_message("threads_job: " + str(threads_job), "yellow")

    ## send for each sample
    print ("+ Calling fastqc for samples...")    
    with concurrent.futures.ThreadPoolExecutor(max_workers=int(max_workers_int)) as executor:
        commandsSent = { executor.submit(fastqc_caller.run_module_fastqc, 
                                         outdir_dict[name[0]], sorted( cluster["sample"].tolist() ), 
                                         name[0], threads_job): name[0] for name, cluster in sample_frame }
        
        for cmd2 in concurrent.futures.as_completed(commandsSent):
            details = commandsSent[cmd2]
            try:
                data = cmd2.result()
            except Exception as exc:
                print ('***ERROR:')
                print (cmd2)
                print('%r generated an exception: %s' % (details, exc))

    print ("+ FASTQC for samples has finished...")    
    
    ## functions.timestamp
    start_time_partial = HCGB_time.timestamp(start_time_partial)

    if (options.skip_report):
        print ("+ No report generation...")
    else:
        print ("\n+ Generating a report using MultiQC module.")
        outdir_report = HCGB_files.create_subfolder("report", outdir)

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
        
        fastqc_report = HCGB_files.create_subfolder("FASTQC", outdir_report)
        fastqc_final_report = HCGB_files.create_subfolder(fold_name, fastqc_report)
        multiQC_report.multiQC_module_call(my_outdir_list, "FASTQC", fastqc_final_report,"")
        print ('\n+ A summary HTML report of each sample is generated in folder: %s' %fastqc_final_report)

    return(outdir_dict)
