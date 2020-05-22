#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez                                      ##
## Copyright (C) 2019-2020 Lauro Sumoy Lab, IGTP, Spain ##
##########################################################
"""
Trimms sequence adapters within fastq reads.
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
from XICRA.scripts import sampleParser
from XICRA.scripts import functions
from XICRA.config import set_config
from XICRA.modules import help_XICRA

##############################################
def run_trimm(options):

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
    
    functions.pipeline_header()
    functions.boxymcboxface("Trimming samples")
    print ("--------- Starting Process ---------")
    functions.print_time()

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
    adapter_dict = {}
    if (options.adapters_a):
        adapters_dict['adapter_a'] = options.adapters_a
    
    if (options.adapters_a):
        adapters_dict['adapter_A'] = options.adapters_A
    
    ## get files
    pd_samples_retrieved = sampleParser.get_files(options, input_dir, "fastq", ("fastq", "fq", "fastq.gz", "fq.gz"))
    
    ## debug message
    if (Debug):
        print (colored("**DEBUG: pd_samples_retrieve **", 'yellow'))
        print (pd_samples_retrieved)

    ## generate output folder, if necessary
    print ("\n+ Create output folder(s):")
    if not options.project:
        functions.create_folder(outdir)
    ## for samples
    outdir_dict = functions.outdir_project(outdir, options.project, pd_samples_retrieved, "trimm")
    
    ## optimize threads
    name_list = set(pd_samples_retrieved["new_name"].tolist())
    threads_job = functions.optimize_threads(options.threads, len(name_list)) ## threads optimization
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
        commandsSent = { executor.submit(cutadapt_caller, sorted(cluster["sample"].tolist()), 
                                         outdir_dict[name], name, threads_job, 
                                         Debug, adapters_dict, options.extra): name for name, cluster in sample_frame }

        for cmd2 in concurrent.futures.as_completed(commandsSent):
            details = commandsSent[cmd2]
            try:
                data = cmd2.result()
            except Exception as exc:
                print ('***ERROR:')
                print (cmd2)
                print('%r generated an exception: %s' % (details, exc))

    print ("\n\n+ Trimming samples has finished...")
    ## functions.timestamp
    start_time_partial = functions.timestamp(start_time_total)

    ## get files generated and generate symbolic link
    if not options.project:
        dir_symlinks = functions.create_subfolder('link_files', outdir)
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
    
        functions.get_symbolic_link(files2symbolic, dir_symlinks)

    if (options.skip_report):
        print ("+ No report generation...")
    else:
        print ("\n+ Generating a report using MultiQC module.")
        outdir_report = functions.create_subfolder("report", outdir)
    
        ## call multiQC report module
        givenList = [ v for v in outdir_dict.values() ]
        my_outdir_list = set(givenList)
        
        ## debug message
        if (Debug):
            print (colored("\n**DEBUG: my_outdir_list for multiqc report **", 'yellow'))
            print (my_outdir_list)
            print ("\n")

        trimm_report = functions.create_subfolder("trimm", outdir_report)
        multiQC_report.multiQC_module_call(my_outdir_list, "Trimmomatic", trimm_report,"")
        print ('\n+ A summary HTML report of each sample is generated in folder: %s' %trimm_report)
        
    print ("\n*************** Finish *******************")
    start_time_partial = functions.timestamp(start_time_total)
    print ("\n+ Exiting trimm module.")
    exit()
    

#############################################
def cutadapt_caller(list_reads, sample_folder, name, threads, Debug, adapters, extra):
    ## check if previously trimmed and succeeded
    filename_stamp = sample_folder + '/.success'
    if os.path.isfile(filename_stamp):
        stamp =    functions.read_time_stamp(filename_stamp)
        print (colored("\tA previous command generated results on: %s [%s -- %s]" %(stamp, name, 'cutadapt'), 'yellow'))
    else:
        # Call cutadapt
        cutadapt_exe = set_config.get_exe('cutadapt')
        code_returned = cutadapt(cutadapt_exe, list_reads, sample_folder, name, threads, Debug, adapters, extra)
        if code_returned:
            functions.print_time_stamp(filename_stamp)
        else:
            print ('** Sample %s failed...' %name)


#############################################
def cutadapt (cutadapt_exe, reads, path, sample_name, num_threads, Debug, adapters, extra):
    """
    
    :param cutadapt_exe:
    :param reads:
    :param path:
    :param sample_name:
    :param num_threads: 
    :param Debug:
    :param adapters
    :param extra:
    
    :type cutadapt_exe:
    :type reads:
    :type path:
    :type sample_name:
    :type num_threads: 
    :type Debug:
    :type adapters: dictionary
    :type extra: string
    
    """
    logfile = os.path.join(path, sample_name + '.cutadapt.log')
    o_param = os.path.join(path, sample_name + '_trim_R1.fastq')
    adapter_3 = ""
    
    if (len(reads) == 2):
        if not adapters['adapter_a'] or adapters_dict['adapter_A']:
             print ("** ERROR: Missing adapter information")
             exit()
        
        p_param = os.path.join(path, sample_name + '_trim_R2.fastq')
        adapter_5 = ""
        ## paired-end mode
        cmd = '%s %s -j %s -m 15 -a %s -A %s -o %s -p %s %s %s > %s' %(cutadapt_exe, extra, 
                                                                       num_threads, adapters['adapter_a'], 
                                                                       adapters_dict['adapter_A'], o_param, 
                                                                       p_param, reads[0], reads[1], logfile)

    elif (len(reads) == 1):
        if not adapters['adapter_a']:
             print ("** ERROR: Missing adapter information")
             exit()
        ## single-end mode:
        cmd = '%s %s -j %s -m 15 -a %s -o %s %s > %s' %(cutadapt_exe, extra, num_threads, 
                                                        adapters['adapter_a'], o_param, reads[0], logfile)    
    else:
        print ('** Wrong number of files provided for sample: %s...' %sample_name)
        return(False)

    ##
    return(functions.system_call(cmd))


    ## cutadapt:
    ## -a adapter_3
    ## -A 3' adapter to be removed from second read in a pair.
    ## 

    