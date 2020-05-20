#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez                                        ##
## Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain        ##
##########################################################
"""
Creates Quality check sequence adapters within fastq reads.
"""
## import useful modules
import os
import sys
import re
import time
from io import open
import shutil
import concurrent.futures
from termcolor import colored
import cutadapt

## import my modules
from XICRA.scripts import multiQC_report
from XICRA.scripts import sampleParser
from XICRA.scripts import functions
from XICRA.config import set_config
from XICRA.modules import help_XICRA

##############################################
def run_QC(options):
    return ()