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
from HCGB import functions
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

    functions.aesthetics_functions.pipeline_header()
    functions.aesthetics_functions.boxymcboxface("Pipeline Configuration")
    print ("--------- Starting Process ---------")
    functions.time_functions.print_time()

    print ("\nCheck dependencies, modules or third party software and print report...")

    ## python version
    functions.aesthetics_functions.print_sepLine("+", 20, False)
    print ('Python:')
    functions.aesthetics_functions.print_sepLine("+", 20, False)

    this_python_version = str(sys.version)
    python_min_version = extern_progs.return_min_version_soft('python')
    if LooseVersion(this_python_version) >= LooseVersion(python_min_version):
        print (colored("Minimum version (%s) satisfied: %s" %( python_min_version, this_python_version), 'green'))
    else:
        print (colored("Minimum version (%s) not satisfied: %s" %(python_min_version, this_python_version), 'red'))
        exit()
        
    ## third-party software
    print ('\n')
    functions.aesthetics_functions.print_sepLine("+", 20, False)
    print ('External dependencies:')
    functions.aesthetics_functions.print_sepLine("+", 20, False)
    
    set_config.check_dependencies(Debug)
    print ('\n')    

    ## python packages
    print ('\n')
    functions.aesthetics_functions.print_sepLine("+", 20, False)
    print ('Python packages:')
    functions.aesthetics_functions.print_sepLine("+", 20, False)

    set_config.check_python_packages(Debug)
    functions.aesthetics_functions.print_sepLine("+", 20, False)
    print ('\n')


    ## R packages
    print ('\n')
    functions.aesthetics_functions.print_sepLine("+", 20, False)
    print ('R packages:')
    functions.aesthetics_functions.print_sepLine("+", 20, False)

    set_config.check_R_packages(Debug)
    functions.aesthetics_functions.print_sepLine("+", 20, False)
    print ('\n')

    