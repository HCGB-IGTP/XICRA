#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez                                        ##
## Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain        ##
##########################################################
"""
Retrieves files within ``other_tools`` directory and returns path to given script specified
"""
## useful imports
import os
from HCGB.functions import files_functions, main_functions

####################################################################
def R_scripts(script, Debug):
    """
    Lists files within ``other_tools/R`` directory and returns path to given script
    """
    RDir = os.path.dirname(os.path.realpath(__file__))
    list_R = main_functions.get_fullpath_list(RDir, Debug)
    
    dict_R = {}
    for f in list_R:
        name = os.path.splitext(os.path.basename(f))[0]
        if (name == script):
            return (f)
        