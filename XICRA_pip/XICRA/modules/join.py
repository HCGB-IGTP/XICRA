#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez                                        ##
## Copyright (C) 2019-2024 Lauro Sumoy Lab, IGTP, Spain        ##
##########################################################
"""
Joins paired-end sequence reads that overlap.
"""
## import useful modules
import os
import time
import concurrent.futures
from termcolor import colored

## import my modules
from XICRA import __version__ as pipeline_version
from XICRA.modules import help_XICRA
from XICRA.config import set_config
from HCGB import sampleParser
import HCGB.functions.aesthetics_functions as HCGB_aes
import HCGB.functions.time_functions as HCGB_time
import HCGB.functions.files_functions as HCGB_files
import HCGB.functions.info_functions as HCGB_info
import HCGB.functions.main_functions as HCGB_main
import HCGB.functions.system_call_functions as HCGB_system

##############################################
def run_join(options):

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
    elif (options.help_join_reads):
        ## information for join reads
        help_XICRA.help_join_reads()
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
    HCGB_aes.boxymcboxface("Join paired-end reads")
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
    
    ## Percentage difference for joining sequences
    if not options.perc_diff:
        options.perc_diff = 0
    
    print ('+ Getting files from input folder... ')
    ## get files
    if options.noTrim:
        print ('+ Mode: fastq.\n+ Extension: ')
        print ("[ fastq, fq, fastq.gz, fq.gz ]\n")
        pd_samples_retrieved = sampleParser.files.get_files(options, input_dir, "fastq", ["fastq", "fq", "fastq.gz", "fq.gz"], options.debug)
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
        HCGB_files.create_folder(outdir)
    ## for samples
    outdir_dict = HCGB_files.outdir_project(outdir, options.project, pd_samples_retrieved, "join", options.debug)
    
    ## optimize threads
    name_list = set(pd_samples_retrieved["new_name"].tolist())
    threads_job = HCGB_main.optimize_threads(options.threads, len(name_list)) ## threads optimization
    max_workers_int = int(options.threads/threads_job)

    ## debug message
    if (Debug):
        print (colored("**DEBUG: options.threads " +  str(options.threads) + " **", 'yellow'))
        print (colored("**DEBUG: max_workers " +  str(max_workers_int) + " **", 'yellow'))
        print (colored("**DEBUG: cpu_here " +  str(threads_job) + " **", 'yellow'))

    print ("+ Joining paired-end sequencing reads for each sample retrieved...")    
    
    # Group dataframe by sample name
    sample_frame = pd_samples_retrieved.groupby(["new_name"])
    
    ## send for each sample
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers_int) as executor:
        commandsSent = { executor.submit(fastqjoin_caller, sorted(cluster["sample"].tolist()), 
                                         outdir_dict[name[0]], name[0], threads_job, options.perc_diff,
                                         Debug): name[0] for name, cluster in sample_frame }

        for cmd2 in concurrent.futures.as_completed(commandsSent):
            details = commandsSent[cmd2]
            try:
                data = cmd2.result()
            except Exception as exc:
                print ('***ERROR:')
                print (cmd2)
                print('%r generated an exception: %s' % (details, exc))

    print ("\n\n+ Joining reads has finished...")
    
    ## TODO: create statistics on joined reads
    ##

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
    runInfo = { "module":"join", "time":time.time(),
                "XICRA version":pipeline_version,
                'sample_info': samples_info,
                'outdir_dict': outdir_dict}
    
    HCGB_info.dump_info_run(info_dir, "join", options, runInfo, options.debug)
    
    ## dump conda details
    HCGB_info.dump_info_conda(info_dir, "join", package_name="XICRA", debug=options.debug)

    print ("\n+ Exiting join module.")

    

    return()

#############################################
def fastqjoin_caller(list_reads, sample_folder, name, threads, perc_diff, Debug):
    ## check if previously joined and succeeded
    filename_stamp = sample_folder + '/.success'
    if os.path.isfile(filename_stamp):
        stamp = HCGB_time.read_time_stamp(filename_stamp)
        print (colored("\tA previous command generated results on: %s [%s -- %s]" %(stamp, name, 'fastqjoin'), 'yellow'))
    else:
        # Call fastqjoin
        fastqjoin_exe = set_config.get_exe('fastqjoin')
        code_returned = fastqjoin(fastqjoin_exe, list_reads, sample_folder, name, threads, perc_diff, Debug)
        if code_returned:
            HCGB_time.print_time_stamp(filename_stamp)
        else:
            print ('** Sample %s failed...' %name)

#############################################
def fastqjoin (fastqjoin_exe, reads, path, sample_name, num_threads, perc_diff, Debug):
    """
    
    :param fastqjoin_exe:
    :param reads:
    :param path:
    :param sample_name:
    :param num_threads: 
    :param Debug:
    
    :type fastqjoin_exe:
    :type reads:
    :type path:
    :type sample_name:
    :type num_threads: 
    :type Debug:
    
    """
    logfile = os.path.join(path, sample_name + '.fastqjoin.log')
    joined_reads = os.path.join(path, sample_name + '_trim_joined.fastq')
    unjoined_1 = os.path.join(path, sample_name + '_trimmed_unjoin_R1.fastq')
    unjoined_2 = os.path.join(path, sample_name + '_trimmed_unjoin_R2.fastq')
    
    ## check paired-end file
    if (len(reads) == 2):
        cmd = fastqjoin_exe + ' -p %s %s %s -o %s -o %s -o %s > %s' %(perc_diff, reads[0], 
                                                                  reads[1], unjoined_1, unjoined_2, 
                                                                  joined_reads, logfile)
    else:
        print ('** Wrong number of files provided for sample: %s...' %sample_name)
        return(False)
    
    

    return(HCGB_system.system_call(cmd))
    