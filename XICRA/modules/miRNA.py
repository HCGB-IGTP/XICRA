#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez                                        ##
## Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain        ##
##########################################################
"""
Create miRNA analysis using sRNAbenchtoolbox and miRTop.
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
from XICRA.scripts import sampleParser
from XICRA.scripts import functions
from XICRA.config import set_config
from XICRA.scripts import parse_gtf
from XICRA.modules import help_XICRA

##############################################
def run_miRNA(options):

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
    
    functions.pipeline_header()
    functions.boxymcboxface("Join paired-end reads")
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
    
    ## get files
    if options.pair:
        pd_samples_retrieved = sampleParser.get_files(options, input_dir, "join", ['_trim_joined'])
    else:
        pd_samples_retrieved = sampleParser.get_files(options, input_dir, "trim", ['_trim'])
    
    ## debug message
    if (Debug):
        print (colored("**DEBUG: pd_samples_retrieve **", 'yellow'))
        print (pd_samples_retrieved)

    ## Additional sRNAbench or miRTop options
    
    ## miRNA_gff: can be set as automatic to download from miRBase
    if not options.miRNA_gff:
        print (colored("** ERROR: No miRNA gff file provided", 'red'))
        exit()
    else:
        options.miRNA_gff = os.path.abspath(options.miRNA_gff)

    ## generate output folder, if necessary
    if not options.project:
        print ("\n+ Create output folder(s):")
        functions.create_folder(outdir)
    ## for samples
    outdir_dict = functions.outdir_project(outdir, options.project, pd_samples_retrieved, "miRNA")
    
    ## optimize threads
    name_list = set(pd_samples_retrieved["new_name"].tolist())
    threads_job = functions.optimize_threads(options.threads, len(name_list)) ## threads optimization
    max_workers_int = int(options.threads/threads_job)

    ## debug message
    if (Debug):
        print (colored("**DEBUG: options.threads " +  str(options.threads) + " **", 'yellow'))
        print (colored("**DEBUG: max_workers " +  str(max_workers_int) + " **", 'yellow'))
        print (colored("**DEBUG: cpu_here " +  str(threads_job) + " **", 'yellow'))

    print ("+ Create a miRNA analysis for each sample retrieved...")    
    
    ## call miRNA_analysis: sRNAbench & miRTop
    # Group dataframe by sample name
    sample_frame = pd_samples_retrieved.groupby(["new_name"])
    
    ## send for each sample
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers_int) as executor:
        commandsSent = { executor.submit(miRNA_analysis, sorted(cluster["sample"].tolist()), 
                                         outdir_dict[name], name, threads_job, options.miRNA_gff,
                                         Debug): name for name, cluster in sample_frame }

        for cmd2 in concurrent.futures.as_completed(commandsSent):
            details = commandsSent[cmd2]
            try:
                data = cmd2.result()
            except Exception as exc:
                print ('***ERROR:')
                print (cmd2)
                print('%r generated an exception: %s' % (details, exc))

    print ("\n\n+ miRNA analysis is finished...")
    
    ## merge all parse gtf files created

    print ("\n*************** Finish *******************")
    start_time_partial = functions.timestamp(start_time_total)
    print ("\n+ Exiting miRNA module.")
    exit()


###############
def miRNA_analysis(reads, folder, name, threads, miRNA_gff, Debug):
    
    ## create sRNAbench
    sRNAbench_folder = functions.create_subfolder('sRNAbench', folder)
    code_success = sRNAbench_caller(reads, sRNAbench_folder, name, threads, Debug) ## Any additional sRNAbench parameter?
        
    if not code_success:
        print ('** miRTop would not be executed for sample %s...' %name)
        return ()
    
    ## create miRTop
    sample_gff = miRTop_caller(sRNAbench_folder, folder, name, threads, miRNA_gff, Debug)
    
    ## parse gtf to accommodate all data
    filename = os.path.join(folder, name + '_XICRA_miRNA.gtf')
    parse_gtf.parse_gtf(sample_gff, filename, name, 'miRNA')

###############       
def sRNAbench_caller(reads, sample_folder, name, threads, Debug):
    # check if previously generated and succeeded
    filename_stamp = sample_folder + '/.success'
    if os.path.isfile(filename_stamp):
        stamp = functions.read_time_stamp(filename_stamp)
        print (colored("\tA previous command generated results on: %s [%s -- %s]" %(stamp, name, 'sRNAbench'), 'yellow'))
    else:
        # Call sRNAbench
        code_returned = sRNAbench(reads, sample_folder, name, threads, Debug)
        if code_returned:
            functions.print_time_stamp(filename_stamp)
        else:
            print ('** Sample %s failed...' %name)
            return(False)
        
    return(True)

###############       
def sRNAbench (reads, outpath, file_name, num_threads, Debug):
    
    sRNAbench_exe = set_config.get_exe("sRNAbench", Debug=Debug)
    sRNAbench_db = os.path.abspath(os.path.join(os.path.dirname(sRNAbench_exe), '..')) ## sRNAtoolboxDB
    logfile = os.path.join(outpath, 'sRNAbench.log')
    
    if (len(reads) > 1):
        print (colored("** ERROR: Only 1 fastq file is allowed please joined reads before...", 'red'))
        exit()
    
    ## create command    
    java_exe = set_config.get_exe('java', Debug=Debug)
    cmd = '%s -jar %s dbPath=%s input=%s output=%s' %(java_exe, sRNAbench_exe, sRNAbench_db, reads[0], outpath)  
    cmd = cmd + ' microRNA=hsa isoMiR=true plotLibs=true graphics=true' 
    cmd = cmd + ' plotMiR=true bedGraphMode=true writeGenomeDist=true'
    cmd = cmd + ' chromosomeLevel=true chrMappingByLength=true > ' + logfile 
    
    return(functions.system_call(cmd))

###############       
def miRTop_caller(sRNAbench_folder, sample_folder, name, threads, miRNA_gff, Debug):
    
    # check if previously generated and succeeded
    mirtop_folder = functions.create_subfolder('miRTop', sample_folder)

    filename_stamp = mirtop_folder + '/.success'
    if os.path.isfile(filename_stamp):
        stamp = functions.read_time_stamp(filename_stamp)
        print (colored("\tA previous command generated results on: %s [%s -- %s]" %(stamp, name, 'miRTop'), 'yellow'))
    else:
        # Call miRTop
        code_returned = miRTop(sRNAbench_folder, mirtop_folder, name, threads, miRNA_gff, Debug)
        if code_returned:
            functions.print_time_stamp(filename_stamp)
        else:
            print ('** Sample %s failed...' %name)
            return(False)
        
        return(True)

###############
def miRTop(sRNAbench_folder, sample_folder, name, threads, miRNA_gff, Debug):

    miRTop_exe = set_config.get_exe('miRTop', Debug=Debug)
    sRNAbench_exe = set_config.get_exe("sRNAbench", Debug=Debug)
    sRNAbench_hairpin = os.path.abspath(os.path.join(os.path.dirname(sRNAbench_exe), '..', 'libs', 'hairpin.fa')) ## sRNAtoolboxDB
    species = 'hsa' #homo sapiens ## set as option if desired
    
    logfile = os.path.join(sample_folder, name + '.log')
    
    ## get sRNAbench info
    reads_annot = os.path.join(sRNAbench_folder, "reads.annotation")
    
    ## check non zero
    if not functions.is_non_zero_file(reads_annot):
        print (colored("\tNo isomiRs detected for sample [%s -- %s]" %(name, 'sRNAbench'), 'yellow'))
        return (False)
    
    ## miRTop analysis
    print ('Creating isomiRs gtf file for sample %s' %name)
    cmd = miRTop_exe + ' gff --sps %s --hairpin %s --gtf %s --format srnabench -o %s %s 2> %s' %(species, sRNAbench_hairpin, miRNA_gff, sample_folder, sRNAbench_folder, logfile)
    
    ## execute
    code_miRTop = functions.system_call(cmd)
    if code_miRTop:
        outdir_stats = functions.create_subfolder('stats', sample_folder)
    
        ## miRTop stats
        print ('Creating isomiRs stats for sample %s' %name)
        cmd_stats = miRTop_exe + ' stats -o %s %s 2>> %s' %(outdir_stats, sample_folder, logfile)
        code_miRTop_stats = functions.system_call(cmd_stats)
    
        ## if both succeeded
        if code_miRTop_stats:
            outdir_gtf = os.path.join(sample_folder, name + '.gff')
            return (outdir_gtf)
    
    ## if any command failed
    return(False)
    