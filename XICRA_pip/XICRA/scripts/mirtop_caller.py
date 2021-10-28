#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez, Marta Lopez & Lauro Sumoy           ##
## Copyright (C) 2019-2021 Lauro Sumoy Lab, IGTP, Spain ##
##########################################################
'''
Calls mirtop to accommodate results
'''
## useful imports
import time
import io
import os
import re
import sys
from sys import argv
from io import open
from termcolor import colored

## import my modules
from HCGB import functions
from XICRA.config import set_config

###############
def miRTop_caller(results_folder, mirtop_folder, name, threads, miRNA_gff, hairpinFasta, format, species, Debug):
    """Checks if the computation has already been performed. If not, it calls
    miRTop()
    :param results_folder: file with the output of the sample from each software 
    :param folder: output miRTop folder
    :param name: sample name
    :param threads: selected threads (by defoult 2)
    :param miRNA_gff: miRNA gff3 annotation file
    :param hairpinFasta: hairpin fasta file
    :param format: 'sRNAbench', 'optimir' or 'seqbuster'
    :param species: species tag ID. Default: hsa (Homo sapiens)
    :param Debug: display complete log.
    :returns: True/False
    """
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
    """Checks if the computation has already been performed. If not, it calls
    miRTop()

    :param results_folder: file with the output of the sample from each software 
    :param folder: output miRTop folder
    :param name: sample name
    :param threads: selected threads (by defoult 2)
    :param miRNA_gff: miRNA gff3 annotation file
    :param hairpinFasta: hairpin fasta file
    :param format: 'sRNAbench', 'optimir' or 'seqbuster'
    :param species: species tag ID. Default: hsa (Homo sapiens)
    :param Debug: display complete log.


    :returns: Folder with the results for the sample in miRTop format
    """

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
    
