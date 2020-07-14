#!/usr/bin/env python3
############################################################
## Jose F. Sanchez                                        ##
## Copyright (C) 2019-2020 Lauro Sumoy Lab, IGTP, Spain   ##
############################################################
"""
Get frequence of reads for each type, variant, etc
"""
## import useful modules
import os
import sys
import re
import time
from io import open
import pandas as pd
from collections import defaultdict
import numpy as np

## import my modules
from XICRA.scripts import functions


## ARGV
if len (sys.argv) < 3:
    print ("\nUsage:")
    print ("python %s freq-table out_file\n" %os.path.realpath(__file__))
    exit()

original_counts = functions.get_data(sys.argv[1], ',', 'index_col=0')
original_counts['total'] = original_counts.sum(axis=1)

print ('##########')
print ('Counts')
print (original_counts)
print ('##########')
print ('')

## get columns
col_list = list(original_counts)
col_list.remove('total')

original_freq = pd.DataFrame()
for miRNA, row in original_counts.iterrows():
    for col in col_list:
        original_freq.loc[miRNA, col] = row[col]/row['total']
    
print ('##########')
print ('Frequence')
print (original_freq)
print ('##########')

df_random = original_freq.apply(np.random.permutation, axis=1, result_type='broadcast')
print (df_random)
print (type(df_random))


