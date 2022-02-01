#!/usr/bin/env python3
############################################################
## Jose F. Sanchez, Marta Lopez & Lauro Sumoy             ##
## Copyright (C) 2019-2021 Lauro Sumoy Lab, IGTP, Spain   ##
############################################################

from HCGB.functions import aesthetics_functions
import HCGB

"""
This module downloads data for genome annotation, miRNA, tRNA and piRNA analysis: 
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
import HCGB.functions.main_functions as HCGB_main
import HCGB.functions.files_functions as HCGB_files
import HCGB.functions.aesthetics_functions as HCGB_aes

from XICRA.config import set_config
from XICRA.modules import help_XICRA
from XICRA.scripts import generate_DE
from XICRA.scripts import MINTMap_caller

##############################################
def run_db(options):
    print()
    
    
##############################################
def miRNA_db(options):
    
    ############################################################
    ## miRNA information: hairpin, mature, str, gff3 from miRBase
    ############################################################
    
    Debug = options.debug
    
    options.miRNA_db = os.path.join(options.database, "miRNA_db")
    functions.files_functions.create_folder(options.miRNA_db)
    
    ## First check if already provided files
    list_files = functions.main_functions.get_fullpath_list(options.miRNA_db, options.debug)
    
    ## Check for files from miRBase:
    miRBase_files = ["hsa.gff3", "hairpin.fa", "mature.fa", "miRNA.str"]
    miRBase_files_dict = {}
    download_data=False
    
    ## get file names
    list_names = list(map(os.path.basename, list_files))
    
    if Debug:
        aesthetics_functions.debug_message("list_files: " + str(list_files))
        aesthetics_functions.debug_message("list_names: " + str(list_names))
    
    for file_req in miRBase_files:
        if Debug:
            aesthetics_functions.debug_message("file_req: " + file_req)
            
        if file_req not in list_names:
            download_data=True
            miRBase_files_dict[file_req] = ""
        else:
            file_retrieved = functions.main_functions.retrieve_matching_files(options.miRNA_db, file_req, options.debug, starts=False)
            if functions.files_functions.is_non_zero_file(file_retrieved[0]):
                miRBase_files_dict[file_req] = file_retrieved[0]
            else:
                miRBase_files_dict[file_req] = ""

    ## If missing, download them, if all files ok, return!
    if Debug:
        aesthetics_functions.debug_message("miRBase_files_dict: ")
        print(miRBase_files_dict)
    
    ## miRBase ftp site
    ftp_site = "https://www.mirbase.org/ftp/CURRENT/"
    
    ## -------------------------------
    ## miRNA_gff: can be set as automatic to download from miRBase
    ## -------------------------------
    if not miRBase_files_dict['hsa.gff3']:
        print ("+ File miRNA gff3 annotation")
        if Debug:
            print (colored("\t** ATTENTION: No miRNA gff file provided", 'yellow'))     
        print (colored("\t** Download it form miRBase", 'green'))
        file_name = options.species + ".gff3"
        ftp_site1 = "https://www.mirbase.org/ftp/CURRENT/genomes/" 
        options.miRNA_gff = functions.main_functions.urllib_request(options.miRNA_db, ftp_site1, file_name, Debug)
        
    else:
        if (options.miRNA_gff):
            print ("\t+ miRNA gff file provided")
            options.miRNA_gff = os.path.abspath(options.miRNA_gff)
        else:
            print ("\t+ miRNA gff file available")
            options.miRNA_gff = miRBase_files_dict['hsa.gff3']

    ## -------------------------------
    ## hairpin: can be set as automatic to download from miRBase
    ## -------------------------------
    if not miRBase_files_dict['hairpin.fa']:
        print ("+ File hairpin fasta")
        if Debug:
            print (colored("\t** ATTENTION: No hairpin fasta file provided", 'yellow'))        
        print (colored("\t** Download it form miRBase", 'green'))
        options.hairpinFasta = functions.main_functions.urllib_request(options.miRNA_db, ftp_site, "hairpin.fa.gz", Debug)
        
    else:
        if (options.hairpinFasta):
            print ("\t+ hairpin fasta file provided")
            options.hairpinFasta = os.path.abspath(options.hairpinFasta)
        else:
            print ("\t+ hairpin fasta file available")
            options.miRNA_gff = miRBase_files_dict['hairpin.fa']
            
   
    ## -------------------------------
    ## mature: can be set as automatic to download from miRBase
    ## -------------------------------
    if not miRBase_files_dict['mature.fa']:
        print ("+ File mature fasta")
        if Debug:
            print (colored("\t** ATTENTION: No mature miRNA fasta file provided", 'yellow'))        
        print (colored("\t** Download it form miRBase", 'green'))
        options.matureFasta = functions.main_functions.urllib_request(options.miRNA_db, ftp_site, "mature.fa.gz", Debug)

    else:
        if (options.matureFasta):
            print ("\t+ mature fasta file provided")
            options.matureFasta = os.path.abspath(options.matureFasta)
        else:
            print ("\t+ miRNA gff file available")
            options.matureFasta = miRBase_files_dict['mature.fa']
        
    ## -------------------------------
    ## miRBase str: can be set as automatic to download from miRBase
    ## -------------------------------
    if not miRBase_files_dict['miRNA.str']:
        print ("+ File miRBase str annotation")
        if Debug:
            print (colored("\t** ATTENTION: No miRBase_str file provided", 'yellow'))        
        print (colored("\t** Download it form miRBase", 'green'))
        options.miRBase_str = functions.main_functions.urllib_request(options.miRNA_db, ftp_site, "miRNA.str.gz", Debug)
        ## extract
        
    else:
        if (options.miRBase_str):
            print ("\t+ miRBase_str file provided")
            options.miRBase_str = os.path.abspath(options.miRBase_str)
        else:
            print ("\t+ miRBase str file available")
            options.miRBase_str = miRBase_files_dict['miRNA.str']
        
    ## -------------------------------
    return (options)    

##############################################
def tRNA_db(database, tRNA_db, debug):
    
    print ("+ Create folder to store several databases: ", database)
    functions.files_functions.create_folder(database)
    
    if not tRNA_db:
        tRNA_db = os.path.join(database, "tRNA_db")
    
    print ("+ Create folder to store tRNA information: ", tRNA_db)
    functions.files_functions.create_folder(tRNA_db)
    
    ## First check if already provided files
    list_files = functions.main_functions.get_fullpath_list(tRNA_db, debug)
    
    ## Check for
    # LookupTable.tRFs.MINTmap_v2.txt
    # OtherAnnotations.MINTmap_v2.txt
    # tRNAspace.Spliced.Sequences.With49ntFlank.MINTmap_v2.fa
    # tables.cfg
    
    ## If missing, download them, if all files ok, return!
    ## TODO    
    
    ## Download information
    ## github repo data: https://github.com/TJU-CMC-Org/MINTmap/tree/release/v2.0-alpha/src/mintmap/mappingbundle/v2
    files = ["https://raw.githubusercontent.com/TJU-CMC-Org/MINTmap/release/v2.0-alpha/src/mintmap/mappingbundle/v2/LookupTable.tRFs.MINTmap_v2.txt",
             "https://raw.githubusercontent.com/TJU-CMC-Org/MINTmap/release/v2.0-alpha/src/mintmap/mappingbundle/v2/OtherAnnotations.MINTmap_v2.txt",
             "https://raw.githubusercontent.com/TJU-CMC-Org/MINTmap/release/v2.0-alpha/src/mintmap/mappingbundle/v2/tRNAspace.Spliced.Sequences.With49ntFlank.MINTmap_v2.fa",
             "https://raw.githubusercontent.com/TJU-CMC-Org/MINTmap/release/v2.0-alpha/src/mintmap/mappingbundle/v2/tables.cfg"
             ]
     
    ## folder generated   
    return (tRNA_db)
    

##############################################
def piRNA_db(database, piRNA_db, debug):
    
    
    print ("+ Create folder to store several databases: ", database)
    functions.files_functions.create_folder(database)
    
    if not piRNA_db:
        piRNA_db = os.path.join(database, "piRNA_db")
    
    print ("+ Create folder to store piRNA information: ", piRNA_db, debug)
    functions.files_functions.create_folder(piRNA_db)
    
    ## First check if already provided files
    list_files = functions.main_functions.get_fullpath_list(piRNA_db, debug)
    
    ## Check for files:
    
    ## piRBase
    # hsa.v3.0.fa
    # hsa.gold.fa
    # hsa.hg19.bed
    
    ## If missing, download them, if all files ok, return!
    ## TODO
    
    return()
    
    ## Download information
    
    ## piRBase
    ## http://bigdata.ibp.ac.cn/piRBase/about.php
    piRBase_hsa = "http://bigdata.ibp.ac.cn/piRBase/download/v3.0/fasta/hsa.v3.0.fa.gz"
    piRBase_hsa_gold = "http://bigdata.ibp.ac.cn/piRBase/download/v3.0/fasta/hsa.gold.fa.gz"
    piRBase_hsa_bed = "http://bigdata.ibp.ac.cn/piRBase/download/v3.0/bed/hsa.hg19.bed.gz"
    
    
    ## some other
    
    ## http://hammelllab.labsites.cshl.edu/software/
    ## TEsmall genome annotation for Homo sapiens hg38
    piRNA="http://labshare.cshl.edu/shares/mhammelllab/www-data/TEsmall/hg38.tar.gz" 
    
    ## DFam database
    ## TSV list of all matches found in the given assembly that score above the GA threshold.
    DFam_db = "https://www.dfam.org/releases/Dfam_3.5/annotations/hg38/hg38.hits.gz"
    
    ## Non-redundant
    ## TSV list of all non-redundant matches found in the given assembly and that score above the GA threshold.
    DFam_db_nrph = "https://www.dfam.org/releases/Dfam_3.5/annotations/hg38/hg38.nrph.hits.gz" 
    
    ## release note DFam: https://www.dfam.org/releases/Dfam_3.5/relnotes.txt
    
    
    return()

####################################################
def piRNA_info(database_folder, species_name="hsa", Debug=False):
    """
    Given a folder it looks for piRNA_db subfolder and gets files
    """
    
    ## to return
    annot_info = {}
    annot_info["piRBase"] = {}
    annot_info["TEsmall_db"] = {}
    
    ## See file in XICRA/config/database/database_structure.txt for tree structure details of this folder
    
    database_folder = os.path.abspath(database_folder)
    piRNA_db = os.path.join( database_folder, "piRNA_db")
    
    ## 
    list_files_piRNA = HCGB_main.get_fullpath_list(piRNA_db, Debug)
    
    ## debugging messages
    if Debug:
        HCGB_aes.debug_message("piRNA_db: " + piRNA_db, "yellow")
        HCGB_aes.debug_message("list_files_piRNA: " + str(list_files_piRNA), "yellow")
        print()
    
    ## loop to retrieve all files required
    for piRNA_files in list_files_piRNA:
        ## debugging messages
        #if Debug:
        #    HCGB_aes.debug_message("piRNA_file: " + piRNA_files, "yellow")
        
        ## get file name
        file_name = os.path.basename(piRNA_files)
        
        ## get gold piRNA for given species
        species_gold_file = str(species_name) + ".gold.fa"
        ## might be some species
        if file_name == species_gold_file:
            if HCGB_files.is_non_zero_file(piRNA_files):
                annot_info["piRBase"]["gold_piRNA"] = piRNA_files
    
        ## get all bed annotation for given species
        if piRNA_files.endswith(".bed"):
            ## store bed files from piRBase
            bed_piRBase = re.search(".*piRBase.*", piRNA_files) ## piRBase db files
            if bed_piRBase:
                if file_name.startswith(species_name) and HCGB_files.is_non_zero_file(piRNA_files):
                        annot_info["piRBase"]["bed_file"] = piRNA_files
                
            ## store bed files from TEsmall_db
            bed_TEsmall_db = re.search(".*TEsmall_db.*", piRNA_files) ## TEsmall_db db files
            if bed_TEsmall_db:
                if species_name=="hsa":
                    
                    ## coding genes
                    if file_name.startswith("exon") and HCGB_files.is_non_zero_file(piRNA_files):
                        annot_info["TEsmall_db"]["exon"] = piRNA_files
                    if file_name.startswith("intron") and HCGB_files.is_non_zero_file(piRNA_files):
                        annot_info["TEsmall_db"]["intron"] = piRNA_files
                    
                    ## miRNA
                    if file_name.startswith("miRNA") and HCGB_files.is_non_zero_file(piRNA_files):
                        annot_info["TEsmall_db"]["miRNA"] = piRNA_files
                    if file_name.startswith("hairpin") and HCGB_files.is_non_zero_file(piRNA_files):
                        annot_info["TEsmall_db"]["hairpin"] = piRNA_files
                    
                    ## TE elements
                    if file_name.startswith("TE") and HCGB_files.is_non_zero_file(piRNA_files):
                        annot_info["TEsmall_db"]["TE"] = piRNA_files
                    
                    ## structural RNA
                    if file_name.startswith("structural_RNA") and HCGB_files.is_non_zero_file(piRNA_files):
                        annot_info["TEsmall_db"]["structural_RNA"] = piRNA_files
        
                    ## piRNA cluster
                    if file_name.startswith("piRNA_cluster") and HCGB_files.is_non_zero_file(piRNA_files):
                        annot_info["TEsmall_db"]["piRNA_cluster"] = piRNA_files

    ## debugging messages
    if Debug:
        HCGB_aes.debug_message("annot_info: " + str(annot_info), "yellow")
        print()
    
    ## concatenate all ncRNA annotation, except piRNA, and coding genes to substract later
    ncRNA_bed = os.path.join(piRNA_db, "ncRNA_merged.bed")
    ncRNA_timestamp = os.path.join(piRNA_db, ".merged_ncRNA_success")
    
    ## check if previously trimmed and succeeded
    if os.path.isfile(ncRNA_timestamp):
        stamp = functions.time_functions.read_time_stamp(ncRNA_timestamp)
        print (colored("\tA previous command generated results on: %s [%s]" %(stamp, 'merged ncRNA'), 'yellow'))
    
    else:
        ## merge
        files_to_merge = [
            annot_info["TEsmall_db"]["exon"], annot_info["TEsmall_db"]["intron"], 
            annot_info["TEsmall_db"]["structural_RNA"], annot_info["TEsmall_db"]["hairpin"], 
            annot_info["TEsmall_db"]["miRNA"]]
        
        ## 
        original_stdout = sys.stdout # Save a reference to the original standard output
        with open(ncRNA_bed, 'w') as f_out:
            sys.stdout = f_out # Change the standard output to the file we created.
            for f_in in files_to_merge:
                for line in open(f_in):
                    line=line.rstrip().split("\t")
                    line[0] = line[0].split("chr")[1]                    
                    print("\t".join(line))

        sys.stdout = original_stdout
        f_out.close()
        HCGB.functions.time_functions.print_time_stamp(ncRNA_timestamp)

    ## save into dictionary info
    annot_info['general'] = {}
    annot_info['general']['ncRNA'] = ncRNA_bed
    annot_info['general']['timestamp'] = ncRNA_timestamp
    
    ## ---------------------------------
    ## return when finished
    return(annot_info)


            
