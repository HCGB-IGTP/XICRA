#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez, Marta Lopez & Lauro Sumoy           ##
## Copyright (C) 2019-2021 Lauro Sumoy Lab, IGTP, Spain ##
##########################################################
'''
Calls MINTMap to create tRNA profile
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
from XICRA.modules import database
import HCGB.functions.aesthetics_functions as HCGB_aes

############################
def download_MINTmap_db(database_folder, tRNA_db, debug):
    """
    Function that calls database module to retrieve tRNA database information required for MINTmap analysis
    """
    ## call database.tRNA_db(database_folder, tRNA_db, debug)
    tRNA_db = database.tRNA_db(database_folder, tRNA_db, debug)
    return (tRNA_db)

############################
def MINTmap_caller(MINTmap_folder, reads, name, num_threads, species, database, Debug):
    # check if previously generated and succeeded
    filename_stamp = MINTmap_folder + '/.success_all'
    if os.path.isfile(filename_stamp):
        stamp = functions.time_functions.read_time_stamp(filename_stamp)
        print (colored("\tA previous command generated results on: %s [%s -- %s]" %(stamp, name, 'MINTmap'), 'yellow'))
        return True

    else:
        # Call MINTMap_analysis
        code_returned = MINTMap_analysis(MINTmap_folder, reads, name, num_threads, species, database, Debug)
        if code_returned:
            functions.time_functions.print_time_stamp(filename_stamp)
            return True
        else:
            print ('** Sample %s failed...' %name)
            return(False)

############################
def MINTMap_analysis(path_folder, reads, name, num_threads, species, database, Debug):
    
    ## check species
    species_code=""
    if species=="hsa":
        species_code="default"
    else:
        print (colored("** ERROR: Not available yet. No mapping bundle available for species " + species, 'red'))
        exit()
    
    ## debug messages 
    if Debug:    
        HCGB_aes.debug_message("species: " + species, "yellow")
        HCGB_aes.debug_message("species_code: " + species_code, "yellow")
    
    ## ATTENTION: MINTmap needs to chdir to output folder
    path_here = os.getcwd()
    
    ## debug messages 
    if Debug:    
        HCGB_aes.debug_message("path_here: " + path_here, "yellow")
            
    filename_stamp = path_folder + '/.success_mintmap'
    if os.path.isfile(filename_stamp):
        stamp = functions.time_functions.read_time_stamp(filename_stamp)
        print (colored("\tA previous command generated results on: %s [%s -- %s]" %(stamp, name, 'MINTmap call'), 'yellow'))
    else:
        # Call MINTMap_analysis
        codeReturn = MINTmap(reads, path_folder, name, num_threads, species_code, database, Debug)
        os.chdir(path_here)
        
        if not codeReturn:
            print ('** Sample %s failed...' %name)
            return False

        ## create time stamp        
        functions.time_functions.print_time_stamp(filename_stamp)

    ## Get MINTmap matrix
    MINTmap_matrix_folder = functions.files_functions.create_subfolder("mintmap_parse", path_folder)

    files = os.listdir(path_folder)
    for item in files:
        abs_path_file = os.path.abspath(os.path.join(path_folder, item))
        
        ## debug messages 
        if Debug:    
            HCGB_aes.debug_message("abs_path_file: " + abs_path_file, "yellow")
        
        ## parse them
        if 'countsmeta' in item:
            continue
        if item.endswith('html'):
            continue
        if 'ambigu' in item:
            amb_file = parse_tRF(abs_path_file, name, MINTmap_matrix_folder, 'amb', Debug)        
        elif 'exclu' in item:
            exc_file = parse_tRF(abs_path_file, name, MINTmap_matrix_folder, 'exc', Debug)
    
    ## debug messages 
    if Debug:    
        HCGB_aes.debug_message("exc_file: " + exc_file, "yellow")
        HCGB_aes.debug_message("amb_file: " + amb_file, "yellow")
        
    
    if functions.files_functions.is_non_zero_file(amb_file) and functions.files_functions.is_non_zero_file(exc_file):
        filename_stamp = path_folder + '/.success_all'
        functions.time_functions.print_time_stamp(filename_stamp)
    
    return(True)

##############
def parse_tRF(pathFile, sample_name, matrix_folder, ident, Debug):
    
    ## tsv file name
    tsv_file = os.path.join(matrix_folder, sample_name + "_" + ident + '.tsv')
    
    ## time stamp
    filename_stamp = matrix_folder + '/.success_parse'
    if os.path.isfile(filename_stamp) and functions.files_functions.is_non_zero_file(tsv_file):
        stamp = functions.time_functions.read_time_stamp(filename_stamp)
        print (colored("\tA previous command generated results on: %s [%s -- %s]" %(stamp, sample_name, 'MINTmap - parse'), 'yellow'))
    else:
        ## Open file
        fil = open(tsv_file, 'w')
        string2write = 'UID\tRead\ttRNA\tvariant\tident\texpression\tsoft\n'
        fil.write(string2write)
        ## Read file
        expression_file = open(pathFile)
        expression_text = expression_file.read()
        expression_lines = expression_text.splitlines()
        
        ## debug messages 
        if Debug:    
            HCGB_aes.debug_message("MINTmap file: " + pathFile, "yellow")
        
        for line in expression_lines:
            # ------------------------------ #
            # Example line:
            # tRF-31-87R8WP9N1EWJ0    TCCCTGGTGGTCTAGTGGTTAGGATTCGGCG    5'-tRF    921    7026.67    452.60    na    trna77_GluCTC_6_+_28949976_28950047@1.31.31, trna80_GluCTC_1_-_161417018_161417089@1.31.31
            # ------------------------------ #
            
            if not line.startswith('#'):
                if not line.startswith('License Plate'):
                    UID = line.split('\t')[0]        ## License Plate
                    seq = line.split('\t')[1]        ## tRF sequence
                    variant = line.split('\t')[2]    ## tRF type
                    expression = line.split('\t')[3] ## unnormalized counts
                    ### there are other RPM counts taking into account several things such as total base pairs, reads, etc. We would use raw counts
                    
                    ## Get tRNA name
                    tRNA_name = line.split('\t')[-1].split(',')[0]
                    tRNA_search = re.search(r"trna.{1,3}\_(.{6})\_(.{1,2})\_.*", tRNA_name)
                    tRNA_family = 'na'
                    if tRNA_search:
                        tRNA_family = tRNA_search.group(1)
                        if (tRNA_search.group(2) == 'MT'):
                            tRNA_family = tRNA_family + '_MT'
                        
                    
                    string2write = UID + '\t' + seq + '\t' + tRNA_family +'\t' + variant +'\t' + ident + '\t' + expression + '\tmintmap\n'
                    
                    ## debug messages 
                    if Debug:    
                        HCGB_aes.debug_message(string2write, "yellow")
                    
                    fil.write(string2write)

        fil.close()

        return(tsv_file)

##############
def MINTmap(reads, outpath, name, num_threads, species_code, database, Debug):
    
    outpath = os.path.abspath(outpath)
    functions.files_functions.create_folder(outpath)
    
    ## change path where the results are required
    os.chdir(outpath)
    
    mintmap_exe = set_config.get_exe("MINTmap", Debug=Debug)
    logfile = os.path.join(outpath, 'MINTmap.log')
    
    ## output
    outpath_file = os.path.join(outpath, name)
    
    if (len(reads) > 1):
        print (colored("** ERROR: Only 1 fastq file is allowed please joined reads before...", 'red'))
        exit()
    
    ## TODO
    ## use -m option with database provided
    
    ## species bundle
    if species_code == "default": 
        ## create command: use default mapping bundle provided with MINTmap 
        cmd = '%s -p %s %s 2> %s' %(mintmap_exe, name, reads[0], logfile)
    else:
        ## create command: use specific mapping bundle path 
        cmd = '%s -p %s -m %s %s 2> %s' %(mintmap_exe, name, species_code, reads[0], logfile)
    
    return(functions.system_call_functions.system_call(cmd))
 