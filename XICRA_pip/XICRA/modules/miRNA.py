#!/usr/bin/env python3
############################################################
## Jose F. Sanchez                                        ##
## Copyright (C) 2019-2020 Lauro Sumoy Lab, IGTP, Spain   ##
############################################################
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
import pandas as pd
from termcolor import colored

## import my modules
from HCGB import sampleParser
from HCGB import functions
from XICRA.config import set_config
from XICRA.modules import help_XICRA
from XICRA.scripts import generate_DE
from XICRA.scripts import reads2tabular

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
    
    functions.aesthetics_functions.pipeline_header()
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
            pd_samples_retrieved = sampleParser.files.get_files(options, input_dir, "fastq", ("fastq", "fq", "fastq.gz", "fq.gz"), options.debug)
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

    ## Additional sRNAbench or miRTop options
    
    ## species
    print ("+ Species provided:", options.species)
    
    ############################################################
    ## miRNA information: hairpin, mature, str, gff3
    ############################################################
    if not (options.database):
        install_path =  os.path.dirname(os.path.realpath(__file__))
        options.database = os.path.join(install_path, "db_files") 
    else:
        options.database = os.path.abspath(options.database)
    
    print ("+ Create folder to store results: ", options.database)
    functions.files_functions.create_folder(options.database)
    
    ## miRNA_gff: can be set as automatic to download from miRBase
    if not options.miRNA_gff:
        print ("+ File miRNA gff3 annotation")
        if Debug:
            print (colored("\t** ATTENTION: No miRNA gff file provided", 'yellow'))     
        print (colored("\t** Download it form miRBase", 'green'))
        file_name = options.species + ".gff3"
        ftp_site = "ftp://mirbase.org/pub/mirbase/CURRENT/genomes/" + file_name 
        options.miRNA_gff = functions.main_functions.urllib_request(options.database, ftp_site, file_name, Debug)
        
    else:
        print ("+ miRNA gff file provided")
        options.miRNA_gff = os.path.abspath(options.miRNA_gff)

    ## hairpin: can be set as automatic to download from miRBase
    if not options.hairpinFasta:
        print ("+ File hairpin fasta")
        if Debug:
            print (colored("\t** ATTENTION: No hairpin fasta file provided", 'yellow'))        
        print (colored("\t** Download it form miRBase", 'green'))
        ftp_site = "ftp://mirbase.org/pub/mirbase/CURRENT/hairpin.fa.gz"
        options.hairpinFasta = functions.main_functions.urllib_request(options.database, ftp_site, "hairpin.fa.gz", Debug)
        
    else:
        print ("+ hairpin fasta file provided")
        options.hairpinFasta = os.path.abspath(options.hairpinFasta)
   
    ## mature: can be set as automatic to download from miRBase
    if not options.matureFasta:
        print ("+ File mature fasta")
        if Debug:
            print (colored("\t** ATTENTION: No mature miRNA fasta file provided", 'yellow'))        
        print (colored("\t** Download it form miRBase", 'green'))
        ftp_site = "ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.gz"
        options.matureFasta = functions.main_functions.urllib_request(options.database, ftp_site, "mature.fa.gz", Debug)

    else:
        print ("+ mature fasta file provided")
        options.matureFasta = os.path.abspath(options.matureFasta)
    
    ## miRBase str: can be set as automatic to download from miRBase
    if not options.miRBase_str:
        print ("+ File miRBase str annotation")
        if Debug:
            print (colored("\t** ATTENTION: No miRBase_str file provided", 'yellow'))        
        print (colored("\t** Download it form miRBase", 'green'))
        ftp_site = "ftp://mirbase.org/pub/mirbase/CURRENT/miRNA.str.gz"
        options.miRBase_str = functions.main_functions.urllib_request(options.database, ftp_site, "miRNA.str.gz", Debug)
        ## extract
        
    else:
        print ("+ miRBase_str file provided")
        options.miRBase_str = os.path.abspath(options.miRBase_str)
    ############################################################
       
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
    
    ## dictionary results
    global results_df
    results_df = pd.DataFrame(columns=("name", "soft", "filename"))
    
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
    
    for soft in soft_list:
        if (soft == "sRNAbench"):
            ## create sRNAbench
            sRNAbench_folder = functions.files_functions.create_subfolder('sRNAbench', folder)
            code_success = sRNAbench_caller(reads, sRNAbench_folder, name, threads, species, Debug) ## Any additional sRNAbench parameter?
                
            if not code_success:
                print ('** miRTop would not be executed for sample %s...' %name)
                return ()
            
            ## create folder for sRNAbench results
            miRTop_folder = functions.files_functions.create_subfolder("sRNAbench_miRTop", folder)
            miRTop_caller(sRNAbench_folder, miRTop_folder, name, threads, miRNA_gff, hairpinFasta, 'sRNAbench', species, Debug)
            
            ## save results in dataframe
            filename = os.path.join(miRTop_folder, 'counts', 'mirtop.tsv')
            results_df.loc[len(results_df)] = name, soft, filename
            
        ###
        if (soft == "optimir"):
            ## create OptimiR analysis
            optimir_folder = functions.files_functions.create_subfolder('OptimiR', folder)
            code_success = optimir_caller(reads, optimir_folder, name, threads, matureFasta, hairpinFasta, miRNA_gff, species, Debug) ## Any additional sRNAbench parameter?
            
            ## create folder for Optimir results
            miRTop_folder = functions.files_functions.create_subfolder("OptimiR_miRTop", folder)
            miRTop_caller(optimir_folder, miRTop_folder, name, threads, miRNA_gff, hairpinFasta, 'optimir', species, Debug)
            
            ## save results in dataframe
            filename = os.path.join(miRTop_folder, 'counts', 'mirtop.tsv')
            results_df.loc[len(results_df)] = name, soft, filename
            
        ###
        if (soft == "miraligner"):
            ## create OptimiR analysis
            miraligner_folder = functions.files_functions.create_subfolder('miraligner', folder)
            code_success = miraligner_caller(reads, miraligner_folder, name, threads, database, species, Debug) 
            
            ## create folder for Optimir results
            miRTop_folder = functions.files_functions.create_subfolder("miraligner_miRTop", folder)
            miRTop_caller(miraligner_folder, miRTop_folder, name, threads, miRNA_gff, hairpinFasta, 'seqbuster', species, Debug)
            
            ## save results in dataframe
            filename = os.path.join(miRTop_folder, 'counts', 'mirtop.tsv')
            results_df.loc[len(results_df)] = name, soft, filename

###############       
def sRNAbench_caller(reads, sample_folder, name, threads, species, Debug):
    # check if previously generated and succeeded
    filename_stamp = sample_folder + '/.success'
    if os.path.isfile(filename_stamp):
        stamp = functions.time_functions.read_time_stamp(filename_stamp)
        print (colored("\tA previous command generated results on: %s [%s -- %s]" %(stamp, name, 'sRNAbench'), 'yellow'))
    else:
        # Call sRNAbench
        code_returned = sRNAbench(reads, sample_folder, name, threads, species, Debug)
        if code_returned:
            functions.time_functions.print_time_stamp(filename_stamp)
        else:
            print ('** Sample %s failed...' %name)
            return(False)
        
    return(True)

###############       
def sRNAbench (reads, outpath, file_name, num_threads, species, Debug):
    
    sRNAbench_exe = set_config.get_exe("sRNAbench", Debug=Debug)
    sRNAbench_db = os.path.abspath(os.path.join(os.path.dirname(sRNAbench_exe), '..')) ## sRNAtoolboxDB
    logfile = os.path.join(outpath, 'sRNAbench.log')
    
    if (len(reads) > 1):
        print (colored("** ERROR: Only 1 fastq file is allowed please joined reads before...", 'red'))
        exit()
    
    ## create command    
    java_exe = set_config.get_exe('java', Debug=Debug)
    cmd = '%s -jar %s dbPath=%s input=%s output=%s' %(java_exe, sRNAbench_exe, sRNAbench_db, reads[0], outpath)  
    cmd = cmd + ' microRNA=%s isoMiR=true plotLibs=true graphics=true' %species 
    cmd = cmd + ' plotMiR=true bedGraphMode=true writeGenomeDist=true'
    cmd = cmd + ' chromosomeLevel=true chrMappingByLength=true > ' + logfile 
    
    return(functions.system_call_functions.system_call(cmd))

###############       
def optimir_caller(reads, sample_folder, name, threads, matureFasta, hairpinFasta, miRNA_gff, species, Debug):
    # check if previously generated and succeeded
    filename_stamp = sample_folder + '/.success'
    if os.path.isfile(filename_stamp):
        stamp = functions.time_functions.read_time_stamp(filename_stamp)
        print (colored("\tA previous command generated results on: %s [%s -- %s]" %(stamp, name, 'OptimiR'), 'yellow'))
    else:
        # Call sRNAbench
	## no species option for OptimiR
        code_returned = optimir(reads, sample_folder, name, threads, matureFasta, hairpinFasta, miRNA_gff,  Debug)
        if code_returned:
            functions.time_functions.print_time_stamp(filename_stamp)
        else:
            print ('** Sample %s failed...' %name)
            return(False)
        
    return(True)

###############       
def optimir (reads, outpath, file_name, num_threads, matureFasta, hairpinFasta, miRNA_gff, Debug):
    
    optimir_exe = set_config.get_exe("optimir", Debug=Debug)
    sRNAbench_db = os.path.abspath(os.path.join(os.path.dirname(optimir_exe), '..')) ## optimir
    logfile = os.path.join(outpath, 'optimir.log')
    errfile = os.path.join(outpath, 'optimir.err')
    
    if (len(reads) > 1):
        print (colored("** ERROR: Only 1 fastq file is allowed please joined reads before...", 'red'))
        exit()
 
    ## create command  
    cmd = "%s process --fq %s --gff_out -o %s --maturesFasta %s --hairpinsFasta %s --gff3 %s > %s 2> %s" %(
        optimir_exe, reads[0], outpath, matureFasta, hairpinFasta, miRNA_gff, logfile, errfile)
    return(functions.system_call_functions.system_call(cmd))

###############       
def miraligner_caller(reads, sample_folder, name, threads, database, species, Debug):
    # check if previously generated and succeeded
    filename_stamp = sample_folder + '/.success'
    if os.path.isfile(filename_stamp):
        stamp = functions.time_functions.read_time_stamp(filename_stamp)
        print (colored("\tA previous command generated results on: %s [%s -- %s]" %(stamp, name, 'sRNAbench'), 'yellow'))
    else:
        # Call miralinger
        code_returned = miraligner(reads, sample_folder, name, database, species, Debug)
        if code_returned:
            functions.time_functions.print_time_stamp(filename_stamp)
        else:
            print ('** Sample %s failed...' %name)
            return(False)
        
    return(True)

###############       
def miraligner (reads, outpath, file_name, database, species, Debug):
    
    miraligner_exe = set_config.get_exe("miraligner", Debug=Debug)
    logfile = os.path.join(outpath, 'miraligner.log')
    
    ## output
    outpath_file = os.path.join(outpath, file_name)
    
    if (len(reads) > 1):
        print (colored("** ERROR: Only 1 fastq file is allowed please joined reads before...", 'red'))
        exit()
    
    ## create tabular information of reads
    tabular_info = os.path.join(outpath, file_name + '-tab.freq.txt')
    reads2tabular.reads2tabular(reads[0], tabular_info)
    
    ## create command 
    java_exe = set_config.get_exe('java', Debug=Debug)
    cmd = '%s -jar %s -db %s -sub 1 -add 3 -trim 3 -s %s -i %s -o %s 2> %s' %(
        java_exe, miraligner_exe, database, species, tabular_info, outpath_file, logfile)
    
    return(functions.system_call_functions.system_call(cmd))


###############       
def miRTop_caller(results_folder, mirtop_folder, name, threads, miRNA_gff, hairpinFasta, format, species, Debug):
    
    # check if previously generated and succeeded
    mirtop_folder_gff = functions.files_functions.create_subfolder('gff', mirtop_folder)
    mirtop_folder_stats = functions.files_functions.create_subfolder('stats', mirtop_folder)
    mirtop_folder_counts = functions.files_functions.create_subfolder('counts', mirtop_folder)
    mirtop_folder_counts = functions.files_functions.create_subfolder('export', mirtop_folder)

    filename_stamp = mirtop_folder_counts + '/.success'
    if os.path.isfile(filename_stamp):
        stamp = functions.time_functions.read_time_stamp(filename_stamp)
        print (colored("\tA previous command generated results on: %s [%s -- %s]" %(stamp, name, 'miRTop'), 'yellow'))
    else:
        # Call miRTop
        code_returned = miRTop(results_folder, mirtop_folder, name, threads, format.lower(), miRNA_gff, hairpinFasta, species, Debug)
        if code_returned:
            functions.time_functions.print_time_stamp(filename_stamp)
        else:
            print ('** Sample %s failed...' %name)
            return(False)
        
        return(True)

###############
def miRTop(results_folder, sample_folder, name, threads, format, miRNA_gff, hairpinFasta, species,Debug):

    miRTop_exe = set_config.get_exe('miRTop', Debug=Debug)
    logfile = os.path.join(sample_folder, name + '.log')
    
    ## folders
    mirtop_folder_gff = os.path.join(sample_folder, 'gff')
    mirtop_folder_stats = os.path.join(sample_folder, 'stats')
    mirtop_folder_counts = os.path.join(sample_folder, 'counts')
    mirtop_folder_export = os.path.join(sample_folder, 'export')
    
    ## get info according to software
    if format == "sRNAbench":
        ## get sRNAbench info
        reads_annot = os.path.join(results_folder, "reads.annotation")
        
        ## check non zero
        if not functions.files_functions.is_non_zero_file(reads_annot):
            print (colored("\tNo isomiRs detected for sample [%s -- %s]" %(name, 'sRNAbench'), 'yellow'))
            return (False)
            
    elif format == "optimir":
        ## get optimir info
        gff3_file = functions.main_functions.retrieve_matching_files(os.path.join(results_folder, "OptimiR_Results"), "gff3", Debug)[0]
        results_folder = gff3_file
        
        ## check non zero
        if not functions.files_functions.is_non_zero_file(gff3_file):
            print (colored("\tNo isomiRs detected for sample [%s -- %s]" %(name, 'optimir'), 'yellow'))
            return (False)
    
    elif format == "seqbuster":
        ## get miraligner info
        mirna_file = functions.main_functions.retrieve_matching_files(results_folder, ".mirna", Debug)[0]
        results_folder = mirna_file
        
        ## check non zero
        if not functions.files_functions.is_non_zero_file(mirna_file):
            print (colored("\tNo isomiRs detected for sample [%s -- %s]" %(name, 'miraligner'), 'yellow'))
            return (False)
    
        
    ## miRTop analysis gff
    filename_stamp_gff = mirtop_folder_gff + '/.success'
    if os.path.isfile(filename_stamp_gff):
        stamp = functions.time_functions.read_time_stamp(filename_stamp_gff)
        print (colored("\tA previous command generated results on: %s [%s -- %s - gff]" %(stamp, name, 'miRTop'), 'yellow'))
    else:
        print ('Creating isomiRs gtf file for sample %s' %name)
        cmd = miRTop_exe + ' gff --sps %s --hairpin %s --gtf %s --format %s -o %s %s 2> %s' %(
                                                    species, hairpinFasta, miRNA_gff, format,
                                                    mirtop_folder_gff, results_folder, logfile)
        
        ## execute
        code_miRTop = functions.system_call_functions.system_call(cmd)
        if code_miRTop:
            functions.time_functions.print_time_stamp(filename_stamp_gff)
        else:
            return(False)
        
    ## miRTop stats
    mirtop_folder_gff_file = os.path.join(mirtop_folder_gff, 'mirtop.gff')

    #filename_stamp_stats = mirtop_folder_stats + '/.success'
    #if os.path.isfile(filename_stamp_stats):
    #    stamp = functions.time_functions.read_time_stamp(filename_stamp_stats)
    #    print (colored("\tA previous command generated results on: %s [%s -- %s - stats]" %(stamp, name, 'miRTop'), 'yellow'))
    #else:
    #    print ('Creating isomiRs stats for sample %s' %name)
    #    cmd_stats = miRTop_exe + ' stats -o %s %s 2>> %s' %(mirtop_folder_stats, mirtop_folder_gff_file, logfile)
    #    code_miRTop_stats = functions.system_call_functions.system_call(cmd_stats)
    #    if code_miRTop_stats:
    #        functions.time_functions.print_time_stamp(filename_stamp_stats)
    #    else:
    #        return(False)
            
    ## miRTop counts
    filename_stamp_counts = mirtop_folder_counts + '/.success'
    if os.path.isfile(filename_stamp_counts):
        stamp = functions.time_functions.read_time_stamp(filename_stamp_counts)
        print (colored("\tA previous command generated results on: %s [%s -- %s - counts]" %(stamp, name, 'miRTop'), 'yellow'))
    else:
        print ('Creating isomiRs counts for sample %s' %name)
        ## if both succeeded
        cmd_stats = miRTop_exe + ' counts -o %s --gff %s --hairpin %s --gtf %s --sps %s 2>> %s' %(mirtop_folder_counts, mirtop_folder_gff_file, hairpinFasta, miRNA_gff, species, logfile)
        code_miRTop_counts = functions.system_call_functions.system_call(cmd_stats)
        
        if code_miRTop_counts:
            functions.time_functions.print_time_stamp(filename_stamp_counts)
        else:
            return(False)
    
    ## miRTop export
    filename_stamp_export = mirtop_folder_export + '/.success'
    if os.path.isfile(filename_stamp_export):
        stamp = functions.time_functions.read_time_stamp(filename_stamp_export)
        print (colored("\tA previous command generated results on: %s [%s -- %s - export]" %(stamp, name, 'miRTop'), 'yellow'))
    else:
        print ('Creating isomiRs export information for sample %s' %name)
        ## if both succeeded
        cmd_export = miRTop_exe + ' export -o %s --hairpin %s --gtf %s --sps %s --format isomir %s 2> %s' %(mirtop_folder_export, hairpinFasta, miRNA_gff, species, mirtop_folder_gff_file, logfile)
        code_miRTop_export = functions.system_call_functions.system_call(cmd_export)
        
        if code_miRTop_export:
            functions.time_functions.print_time_stamp(filename_stamp_export)
        else:
            return(False)
    
    ## return all success
    outdir_tsv = os.path.join(mirtop_folder_counts, 'mirtop.tsv')
    return (outdir_tsv)
    
    ## if any command failed
    #return(False)
    
