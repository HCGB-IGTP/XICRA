#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez                                        ##
## Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain        ##
##########################################################
"""
Calls multiQC to generate HTML statistics reports.
"""
## useful imports
import os
import io
import sys
from io import open
from sys import argv
from termcolor import colored

## import my modules
from HCGB import functions
from XICRA.config import set_config

############
def multiQC_module_call(givenList, name, path, option):
    """
    Prepares files for multiQC report generation.
    
    :param givenList: List of folder to search for multiQC report.
    :param name: Name to include in the html report.
    :param path: Absolute path for the output folder.
    :param option: Some options to provide to multiQC_call.
    
    :type givenList: list
    :type name: string
    :type path: string
    :type option: string
    
    .. seealso:: This function depends on other XICRA functions called:
    
        - :func:`XICRA.scripts.functions.printList2file`
        
        - :func:`XICRA.scripts.multiQC_report.multiQC_call`
    
    """
    pathFile = path + '/' + 'samples.txt'
    functions.main_functions.printList2file(pathFile, givenList)
    multiQC_call(pathFile, name, path, option)    
    
############
def multiQC_call(pathFile, name, folder, option):
    """
    multiQC_ report generation call.
    
    :param pathFile: File containing list of files to include in report.
    :param name: Name to include in the html report.
    :param folder: Absolute path for the output folder.
    :param option: Options to provide to multiQC call.
    
    :type pathFile: string
    :type name: string 
    :type folder: string 
    :type option: string
    
    :returns: :func:`XICRA.scripts.functions.system_call_functions.system_call` output (OK/FALSE)
        
    .. seealso:: This function depends on other XICRA functions called:
    
        - :func:`XICRA.scripts.functions.system_call_functions.system_call`
    
    """
    multiqc_bin = set_config.get_exe("multiqc")
    ## set options for call
    cmd = "%s --force -o %s -n %s -l %s -p -i 'MultiQC report' -b 'HTML report generated for multiple samples and steps' %s" %(multiqc_bin, folder, name, pathFile, option)
    
    ## if a report was previously generated in the folder 
    ## force to delete and generate a new one
    return(functions.system_call_functions.system_call(cmd))
