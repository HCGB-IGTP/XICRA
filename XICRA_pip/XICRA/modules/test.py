#!/usr/bin/env python3

"""
test example file
"""

import os
import HCGB.functions
import shutil

## Download from git subset2test folder and run sh command

def run_test(options):
    print ("** Downloading folder containing test reads from github")
    
    ## Download repo
    str_command = "gitdir https://github.com/HCGB-IGTP/XICRA/tree/master/subset2test/"
    HCGB.functions.system_call_functions.system_call(str_command)

    shutil.move("./subset2test/subset_PE/", "./subset_PE/")
    shutil.move("./subset2test/subset_SE/", "./subset_SE/")
    shutil.move("./subset2test/subset_tRNA/", "./subset_tRNA/")
    shutil.move("./subset2test/test_subset.sh", "./test_subset.sh")
    
    shutil.rmtree("./subset2test/")
    
    print ("** Now, either execute sh test_subset.sh or step by step")
    
