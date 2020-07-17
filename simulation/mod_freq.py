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
import random

## import my modules
from XICRA.scripts import functions

## get frequencies
def get_freq(given_df, col_list):
    df_freq = pd.DataFrame()
    for miRNA, row in given_df.iterrows():
        for col in col_list:
            df_freq.loc[miRNA, col] = row[col]/row['total']
    return (df_freq)

## ARGV
if len (sys.argv) < 3:
    print ("\nUsage:")
    print ("python %s freq-table out_file debug_bool\n" %os.path.realpath(__file__))
    exit()

## original counts
original_counts = functions.get_data(sys.argv[1], ',', 'index_col=0')
col_list = list(original_counts) ## get columns

## add missing data
modified_counts = original_counts.copy(deep=True)
for col in col_list:
    modified_counts.loc[modified_counts.sample(frac=0.35).index, col] = pd.np.nan

## randomize 
random_counts = modified_counts.apply(np.random.permutation, axis=1, result_type='broadcast')
random_counts[np.isnan(random_counts)] = 0
random_counts['total'] = random_counts.sum(axis=1)

## get frequencies
random_freqs = get_freq(random_counts, col_list)

if (sys.argv[3]):
    print ('##########')
    print ('Counts')
    print (original_counts)
    print ('##########')
    print ('')

    print ('##########')
    print ('Frequence')
    original_counts['total'] = original_counts.sum(axis=1)
    original_freq = get_freq(original_counts, col_list)
    print (original_freq)
    print ('##########')

    ## print randomize counts & frequencies
    print ('##########')
    print ('Random Counts')
    print (random_counts)
    print ('##########')
    print ('')

    print ('##########')
    print ('Frequence')
    print (random_freqs)
    print ('##########')

## adjust to 100
new_random = pd.DataFrame(columns=col_list)
for miRNA, row in random_freqs.iterrows():
    for col in col_list:
        if row[col]==0:
            new_random.loc[miRNA, col] = 0
        else:
            new_random.loc[miRNA, col] = int(row[col]*100) 

new_random['total'] = new_random.sum(axis=1)
for miRNA, row in new_random.iterrows():
    if row['total']!=100:
        sum = 100 - int(row['total']) 
        rnd = random.sample(col_list, 1)
        row[rnd] += sum

new_random = new_random.drop(columns=['total'])
new_random['total'] = new_random.sum(axis=1)


print ('##########')
print ('Counts')
print (original_counts)
print ('##########')
print ('')

## print randomize counts & frequencies
print ('##########')
print ('Counts adjusted')
print (new_random)
print ('##########')
print ('')