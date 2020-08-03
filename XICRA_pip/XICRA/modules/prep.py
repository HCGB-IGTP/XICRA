#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez                                        ##
## Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain        ##
##########################################################
"""
Get fastq files and prepare them for other purposes. 
"""

## import useful modules
import time
import os
import sys
from io import open
import shutil
import pandas as pd
from termcolor import colored

## import my modules
from XICRA.modules import help_XICRA
from HCGB import functions
from HCGB import sampleParser

################################
def run_prep(options):
    """
    Main function of the prep module.
    
    This module prepares fastq files for later usage. It initially checks the length
    of the name and advises the user to rename samples if exceeded. 
    
    This module allows to user to copy files into the project folder initiate or only link using
    a symbolic link to avoid duplicated raw data. 

    """
    
    ## help_format option
    if (options.help_format):
        help_XICRA.help_fastq_format()
        exit()
        
    functions.aesthetics_functions.pipeline_header()
    functions.aesthetics_functions.boxymcboxface("Preparing samples")
    print ("--------- Starting Process ---------")
    functions.time_functions.print_time()
    
    ## init time
    start_time_total = time.time()
    
    ## absolute path for in & out
    input_dir = os.path.abspath(options.input)
    outdir = os.path.abspath(options.output_folder)

    ### set as default paired_end mode
    if (options.single_end):
        options.pair = False
    else:
        options.pair = True

    ## output folder    
    print ("\n+ Create output folder(s):")
    functions.files_functions.create_folder(outdir)

    ## project
    if options.detached:
        options.project = False
    else:
        options.project = True
    
    ## default options
    if not options.include_lane:
        options.include_lane = False
    
    if not options.include_all:
        options.include_all = False

    ### info
    final_dir = ""
    if (options.project):
        print ("+ Generate a directory containing information within the project folder provided")
        final_dir = functions.files_functions.create_subfolder("info", outdir)
    else:
        final_dir = outdir
    
    ## get files
    print ()
    functions.aesthetics_functions.print_sepLine("-",50, False)
    print ('+ Getting files from input folder... ')
    print ('+ Mode: fastq.\n+ Extension: ')
    print ("[ fastq, fq, fastq.gz, fq.gz ]\n")
    
    pd_samples_retrieved = sampleParser.files.get_files(options, input_dir, "fastq", ("fastq", "fq", "fastq.gz", "fq.gz"), options.debug)
        
    ## Information returned in pd_samples_retrieved
    ### sample, dirname, name, name_len, lane, read_pair, lane_file, ext, gz
    
    if options.debug:
        print (colored("** DEBUG: pd_samples_retrieved", 'yellow'))
        print (pd_samples_retrieved)
        #functions.print_all_pandaDF(pd_samples_retrieved)
    
    ## time stamp
    start_time_partial = functions.time_functions.timestamp(start_time_total)
    
    ## check character limitation
    list_lengths = pd_samples_retrieved.loc[:,'name_len'].to_list()
    if any(i > 25 for i in list_lengths):
        print (colored("\t ** Name lengths exceeds the 25 character limitation...", 'yellow'))
        if not (options.rename):
            print (colored("** ERROR: Rename files or provide --rename option...", 'red'))
            exit()

    ### rename files 
    if (options.rename):
        options.rename = os.path.abspath(options.rename)
        if not functions.files_functions.is_non_zero_file(options.rename):
            print (colored("** ERROR: File provided with rename information is not readable.", 'red'))
            print (options.rename)
            exit()
            
        names_retrieved = pd.read_csv(options.rename, sep=',', 
                                    index_col=0, squeeze=True, 
                                    header=None).to_dict() ## read csv to dictionary
        if (options.debug):
            print (colored('** DEBUG: names_retrieved', 'yellow'))
            print (names_retrieved)
            
        ## TODO: check integrity of new names and special characters
    
        ## print to a file
        timestamp = functions.time_functions.create_human_timestamp()
        rename_details = final_dir + '/' + timestamp + '_prep_renameDetails.txt'
        rename_details_hd = open(rename_details, 'w')
    
        ## rename files         
        for index, row in pd_samples_retrieved.iterrows():
            if (row['gz']):
                extension_string = row['ext'] + row['gz']
            else:
                extension_string = row['ext']
            
            if options.single_end:
                renamed = names_retrieved[row['name']] + '.' + extension_string
            else:
                renamed = names_retrieved[row['name']] + '_' + row['read_pair'] + '.' + extension_string
            
            ## modify frame
            pd_samples_retrieved.loc[index, 'new_file'] = renamed
            pd_samples_retrieved.loc[index, 'new_name'] = names_retrieved[row['name']]
            ## save in file
            string = row['sample'] + '\t' + renamed + '\n'
            rename_details_hd.write(string)
            
            if (options.debug):
                print (colored('** DEBUG: rename', 'yellow'))
                print ("Original: ", row['name'])
                print ("Renamed: ", names_retrieved[row['name']])
                print ("File:", renamed)
        
        rename_details_hd.close()    

        ##elif (options.single_end): It should work for both
        print ("+ Sample files have been renamed...")
    else:
        pd_samples_retrieved['new_file'] = pd_samples_retrieved['file']
 
    ## create outdir for each sample
    outdir_dict = functions.files_functions.outdir_project(outdir, options.project, pd_samples_retrieved, "raw", options.debug)    
        
    ## merge option
    if (options.merge_Reads):
        print ("+ Sample files will be merged...")
        ## TODO: check when rename option provided
        pd_samples_merged = sampleParser.merge.one_file_per_sample(
            pd_samples_retrieved, outdir_dict, options.threads,    
            final_dir, options.debug)
        
        if (options.rename):
            print ("+ Merge files have been renamed...")
        else:
            print ("+ Sample files have been merged...")
        
        ## process is finished here
        print ("\n*************** Finish *******************")
        start_time_partial = functions.time_functions.timestamp(start_time_total)
    
        print ("+ Exiting prep module.")
        exit()
    
    ## debugging messages
    if (options.debug):
        print (colored("** DEBUG: pd_samples_retrieved", 'yellow'))
        #functions.print_all_pandaDF(pd_samples_retrieved)
        print (pd_samples_retrieved)
        print (colored("** DEBUG: outdir_dict", 'yellow'))
        print (outdir_dict)
    
    ## copy or create symbolic link for files
    if (options.copy_reads):
        print ("+ Sample files will be copied...")
        ## print to a file
        timestamp = functions.time_functions.create_human_timestamp()
        copy_details = final_dir + '/' + timestamp + '_prep_copyDetails.txt'
        copy_details_hd = open(copy_details, 'w')
    else:
        print ("+ Sample files will be linked...")    
    
    list_reads = []
    for index, row in pd_samples_retrieved.iterrows():
        if (options.copy_reads):
            ## TODO: debug & set threads to copy faster
            shutil.copy(row['sample'], os.path.join(outdir_dict[row['new_name']], row['new_file'] ))            
            string = row['sample'] + '\t' + os.path.join(outdir_dict[row['new_name']], row['new_file']) + '\n'
            copy_details_hd.write(string)            
        else:
            list_reads.append(row['new_file'])
            
            if options.project:
                functions.files_functions.get_symbolic_link_file(row['sample'], 
                                                 os.path.join(outdir_dict[row['new_name']], row['new_file']))

    if (options.copy_reads):
        print ("+ Sample files have been copied...")
        copy_details_hd.close()
    else:
        if not options.project:
            functions.files_functions.get_symbolic_link(list_reads, outdir)
    
    print ("\n*************** Finish *******************")
    start_time_partial = functions.time_functions.timestamp(start_time_total)

    print ("+ Exiting prep module.")
    return()
