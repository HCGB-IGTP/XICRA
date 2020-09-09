#!/usr/bin/env python3
############################################################
## Jose F. Sanchez                                        ##
## Copyright (C) 2019-2020 Lauro Sumoy Lab, IGTP, Spain   ##
############################################################
## useful imports
import time
import io
import os
import re
import sys
from sys import argv
import subprocess

from HCGB.functions import system_call_functions
from HCGB.functions import files_functions

############################################################
def load_Genome(folder, STAR_exe, genomeDir, num_threads):
    
    ## --genomeLoad LoadAndExit
    LoadDir = 'LoadMem'
    Load_folder = files_functions.create_subfolder(LoadDir, folder)
    cmd_LD = "%s --genomeDir %s --runThreadN %s --outFileNamePrefix %s --genomeLoad LoadAndExit" %(
        STAR_exe, genomeDir, num_threads, Load_folder)
    
    print ('\t+ Loading memory for STAR mapping')
    load_code = system_call_functions.system_call(cmd_LD, False, True)
    return (load_code)

############################################################
def remove_Genome(STAR_exe, genomeDir, remove_folder, num_threads):
    
    ## --genomeLoad Remove
    removeDir = 'RemoveMem'
    remove_folder = files_functions.create_subfolder(removeDir, folder)
    cmd_RM = "%s --genomeDir %s --outFileNamePrefix %s --runThreadN %s --genomeLoad Remove" %(
        STAR_exe, genomeDir, folder, num_threads)
    
    ## send command    
    print ('\t+ Removing memory loaded for STAR mapping')
    remove_code = system_call_functions.system_call(cmd_RM, False, True)
    return (remove_code)

############################################################
def mapReads(option, read, folder, name, STAR_exe, genomeDir, limitRAM_option, num_threads):
    ## all this options and parameters have been obtained 
    ## from https://www.encodeproject.org/rna-seq/small-rnas/

    ## open file
    print("\t + Mapping sample %s using STAR" %name)
    
    out_folder = files_functions.create_subfolder(name, folder)
    logfile = os.path.join(out_folder, 'STAR.log')
    errfile = os.path.join(out_folder, 'STAR.err')
    bam_file = os.path.join(out_folder, 'Aligned.sortedByCoord.out.bam')

    ## prepare command
    cmd = "%s --genomeDir %s --runThreadN %s --readFilesIn %s " %(STAR_exe, genomeDir, num_threads, jread)
    cmd = cmd + "--limitBAMsortRAM %s --outFileNamePrefix %s" %(limitRAM_option, sample_folder)

    ## some common options
    cmd = cmd + "--alignSJDBoverhangMin 1000 --outFilterMultimapNmax 1 --outFilterMismatchNoverLmax 0.03 "
    cmd = cmd + "--outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 16 "
    cmd = cmd + "--alignIntronMax 1 --outSAMheaderHD @HD VN:1.4 SO:coordinate --outSAMtype BAM SortedByCoordinate "
    
    ## Multiple samples or just one?
    if option == 'LoadAndKeep':
        cmd = cmd + "--genomeLoad LoadAndKeep "
    else:
        cmd = cmd ## ???

    ## logfile & errfile
    cmd = cmd + ' > ' + logfile + ' 2> ' + errfile
    
    ## sent command
    mapping_code = system_call_functions.system_call(cmd, False, True)
    return mapping_code

###############

###########
def main():
    
    
    return 

######
if __name__== "__main__":
    main()