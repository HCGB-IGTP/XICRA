#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez                                      ##
## Copyright (C) 2019-2020 Lauro Sumoy Lab, IGTP, Spain ##
##########################################################
"""
Trims sequence adapters within fastq reads.
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
from XICRA.scripts import multiQC_report
from XICRA.scripts import cutadapt_caller
from XICRA.config import set_config
from XICRA.modules import help_XICRA
from XICRA.modules import qc
from HCGB import functions
from HCGB import sampleParser

##############################################
def run_trim(options):
    """Main function of the module, organizes the trimming process.

    First, checks if the adapter(s) sequence(s) have been provided by the user:
    - adapters_a or adapters_A
    If there is no adapter sequence provided the process will be stopped.  

    If the adapters have been introduced, it calls cutadapt_caller() for each sample in parallel.
    Finally, generates a report using MultiQC module if desired.
    
    :param options: input parameters introduced by the user. See XICRA trim -h.

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
    elif (options.help_trimm_adapters):
        ## help on trimm adapters
        help_XICRA.print_help_adapters()
        exit()
    elif (options.help_project):
        ## information for project
        help_XICRA.project_help()
        exit()
    elif (options.help_multiqc):
        ## information for Multiqc
        help_XICRA.multiqc_help()
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
    functions.aesthetics_functions.boxymcboxface("Trimming samples")
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
    
    # Trimming adapters

    ## check adapters provided
        ## options.adapters_a
        ## options.adapters_A
        ## options.extra
        
    ## no adapters provided
    if (not options.adapters_a and not options.adapters_A and not options.extra):
        print (colored("** ERROR: No adapter trimming options provided...", 'red'))
        print ("Please provide any option")
        exit()
    
    ## create dictionary with 
    adapters_dict = {}
    if (options.adapters_a):
        adapters_dict['adapter_a'] = options.adapters_a
    
    if (options.adapters_a):
        adapters_dict['adapter_A'] = options.adapters_A
    
    ## set default
    #if not options.min_len_read:
    #    options.min_len_read=15
    
    ## get files
    print ('+ Getting files from input folder... ')
    print ('+ Mode: fastq.\n+ Extension: ')
    print ("[ fastq, fq, fastq.gz, fq.gz ]\n")
    pd_samples_retrieved = sampleParser.files.get_files(options, input_dir, "fastq", ["fastq", "fq", "fastq.gz", "fq.gz"], options.debug)
    
    ## debug message
    if (Debug):
        print (colored("**DEBUG: pd_samples_retrieve **", 'yellow'))
        print (pd_samples_retrieved)

        print (colored("**DEBUG: adapters_dict **",'yellow'))
        print (adapters_dict)

    ## generate output folder, if necessary
    print ("\n+ Create output folder(s):")
    if not options.project:
        functions.files_functions.create_folder(outdir)
    ## for samples
    outdir_dict = functions.files_functions.outdir_project(outdir, options.project, pd_samples_retrieved, "trimm", options.debug)
    
    ## optimize threads
    name_list = set(pd_samples_retrieved["new_name"].tolist())
    threads_job = functions.main_functions.optimize_threads(options.threads, len(name_list)) ## threads optimization
    max_workers_int = int(options.threads/threads_job)

    ## debug message
    if (Debug):
        print (colored("**DEBUG: options.threads " +  str(options.threads) + " **", 'yellow'))
        print (colored("**DEBUG: max_workers " +  str(max_workers_int) + " **", 'yellow'))
        print (colored("**DEBUG: cpu_here " +  str(threads_job) + " **", 'yellow'))

    print ("+ Trimming adapters for each sample retrieved...")    
    
    # Group dataframe by sample name
    sample_frame = pd_samples_retrieved.groupby(["new_name"])
    
    ## send for each sample
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers_int) as executor:
        commandsSent = { executor.submit(cutadapt_caller.caller, sorted(cluster["sample"].tolist()), 
                                         outdir_dict[name], name, threads_job, 
                                         options.min_read_len, Debug, adapters_dict, options.extra): name for name, cluster in sample_frame }

        for cmd2 in concurrent.futures.as_completed(commandsSent):
            details = commandsSent[cmd2]
            try:
                data = cmd2.result()
            except Exception as exc:
                print ('***ERROR:')
                print (cmd2)
                print('%r generated an exception: %s' % (details, exc))

    print ("\n\n+ Trimming samples has finished...")
    
    ## functions.time_functions.timestamp
    start_time_partial = functions.time_functions.timestamp(start_time_total)

    ## get files generated and generate symbolic link
    if not options.project:
        dir_symlinks = functions.files_functions.create_subfolder('link_files', outdir)
        files2symbolic = []
        folders = os.listdir(outdir)

        ## debug message
        if (Debug):
            print (colored("**DEBUG: generate symbolic links for each file in " + dir_symlinks + "**", 'yellow'))
        
        for fold in folders:
            if fold.endswith(".log"):
                continue
            else:
                this_folder = outdir + '/' + fold
                subfiles = os.listdir(this_folder)
                for files in subfiles:
                    files_search = re.search(r".*trim_R\d{1}.*", files) ## only paired-end. Todo: single end
                    if files_search:
                        files2symbolic.append(this_folder + '/' + files)
    
        functions.files_functions.get_symbolic_link(files2symbolic, dir_symlinks)

    if (options.skip_report):
        print ("+ No report generation...")
    else:
        print ("\n+ Generating a report using MultiQC module.")
        outdir_report = functions.files_functions.create_subfolder("report", outdir)
    
        ## call multiQC report module
        givenList = [ v for v in outdir_dict.values() ]
        my_outdir_list = set(givenList)
        
        ## debug message
        if (Debug):
            print (colored("\n**DEBUG: my_outdir_list for multiqc report **", 'yellow'))
            print (my_outdir_list)
            print ("\n")

        trimm_report = functions.files_functions.create_subfolder("trim", outdir_report)
        multiQC_report.multiQC_module_call(my_outdir_list, "Cutadapt", trimm_report,"")
        print ('\n+ A summary HTML report of each sample is generated in folder: %s' %trimm_report)
        
        ## QC analysis for trimmed reads
        if (Debug):
            print (colored("** Beginning FAStQC analysis **", 'red'))

        ## functions.time_functions.timestamp
        start_time_partial = functions.time_functions.timestamp(start_time_partial)

    ## create FASTQC calling for trimmed reads
    pd_samples_retrieved_trimmed = sampleParser.files.get_files(options, input_dir, "trim", ['_trim'], options.debug)
    qc.fastqc(pd_samples_retrieved_trimmed, outdir, options, "trimmed", start_time_partial, Debug)
        
    print ("\n*************** Finish *******************")
    start_time_partial = functions.time_functions.timestamp(start_time_total)
    print ("\n+ Exiting trim module.")
    exit()

