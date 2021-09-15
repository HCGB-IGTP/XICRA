#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez                                      ##
## Copyright (C) 2019-2020 Lauro Sumoy Lab, IGTP, Spain ##
##########################################################

## useful imports
import time
import io
import os
import sys
from termcolor import colored
from distutils.version import LooseVersion

## import my modules
import HCGB.functions.time_functions as HCGB_time
import HCGB.functions.aesthetics_functions as HCGB_aes
from XICRA.config import extern_progs
from XICRA.config import set_config

################
def run_config(options):
    ## init time
    start_time_total = time.time()

    ## debugging messages
    global Debug
    if (options.debug):
        Debug = True
    else:
        Debug = False

    HCGB_aes.pipeline_header('XICRA')
    HCGB_aes.boxymcboxface("Pipeline Configuration")
    print ("--------- Starting Process ---------")
    HCGB_time.print_time()

    print ("\nCheck dependencies, modules or third party software and print report...")

    ## python version
    HCGB_aes.print_sepLine("+", 20, False)
    print ('Python:')
    HCGB_aes.print_sepLine("+", 20, False)

    this_python_version = str(sys.version)
    python_min_version = extern_progs.return_min_version_soft('python')
    if LooseVersion(this_python_version) >= LooseVersion(python_min_version):
        print (colored("Minimum version (%s) satisfied: %s" %( python_min_version, this_python_version), 'green'))
    else:
        print (colored("Minimum version (%s) not satisfied: %s" %(python_min_version, this_python_version), 'red'))
        exit()
        
    ## third-party software
    print ('\n')
    HCGB_aes.print_sepLine("+", 20, False)
    print ('External dependencies:')
    HCGB_aes.print_sepLine("+", 20, False)
    
    set_config.check_dependencies(Debug)
    print ('\n')    

    ## python packages
    print ('\n')
    HCGB_aes.print_sepLine("+", 20, False)
    print ('Python packages:')
    HCGB_aes.print_sepLine("+", 20, False)

    set_config.check_python_packages(Debug)
    HCGB_aes.print_sepLine("+", 20, False)
    print ('\n')


    ## R packages
    print ('\n')
    HCGB_aes.print_sepLine("+", 20, False)
    print ('R packages:')
    HCGB_aes.print_sepLine("+", 20, False)

    set_config.check_R_packages(Debug)
    HCGB_aes.print_sepLine("+", 20, False)
    print ('\n')

    