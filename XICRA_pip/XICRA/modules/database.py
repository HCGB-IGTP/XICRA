#!/usr/bin/env python3
############################################################
## Jose F. Sanchez, Marta Lopez & Lauro Sumoy             ##
## Copyright (C) 2019-2021 Lauro Sumoy Lab, IGTP, Spain   ##
############################################################
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
    
    ## First check if already provided files
    list_files = functions.main_functions.get_fullpath_list(options.database, options.debug)
    
    ## Check for
    # hsa.gff3
    # hairpin.fa
    # mature.fa
    # miRNA.str
    
    ## If missing, download them, if all files ok, return!
    ## TODO
    
    ## -------------------------------
    ## miRNA_gff: can be set as automatic to download from miRBase
    ## -------------------------------
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

    ## -------------------------------
    ## hairpin: can be set as automatic to download from miRBase
    ## -------------------------------
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
   
    ## -------------------------------
    ## mature: can be set as automatic to download from miRBase
    ## -------------------------------
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
    
    ## -------------------------------
    ## miRBase str: can be set as automatic to download from miRBase
    ## -------------------------------
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
    
    print ("+ Create folder to store tRNA information: ", piRNA_db, debug)
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
    
