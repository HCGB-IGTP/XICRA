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
import argparse

## import my modules
from XICRA.scripts import functions
from XICRA.scripts import reads2tabular

## get frequencies
def get_freq(given_df, col_list):
    df_freq = pd.DataFrame()
    for miRNA, row in given_df.iterrows():
        for col in col_list:
            if row[col]==0:
                df_freq.loc[miRNA, col] = 0
            else:
                df_freq.loc[miRNA, col] = row[col]/row['total']
    return (df_freq)


#####################################################
parser = argparse.ArgumentParser(prog='mod_freq.py', formatter_class=argparse.RawDescriptionHelpFormatter, 
                                 description='''

mod_freq.py: Modified given frequencies and select isomiRs

Version: 0.1
License: GPLv3

USAGE: python mod_freq.py --freq table.freq.csv --out out_name [--debug] 
''', epilog="Original code: JFSanchezHerrero")

#####################################################
parser.add_argument('-f', '--freq', action='store', help='Table with original variant frequencies to modify', required=True)
parser.add_argument('-o', '--out', action='store', help='Output names', required=True)
parser.add_argument('--debug', action='store_true', default=False, help='Developer messages')
parser.add_argument('--random_rows', action='store', type=int, help='Numbers of miRNA to subset', default=100)
args = parser.parse_args()
#####################################################

## original counts
print ("# Read original frequency table")
original_counts = functions.get_data(args.freq, ',', 'index_col=0')
col_list = list(original_counts) ## get columns

## drop NAs
print ("# Remove any rows containing NAs from frequency table")
original_counts = original_counts.dropna()

## subset 100 rows
print ("# Randomly subsetting rows")
subset_df = original_counts.sample(n=args.random_rows)

## add missing data
print ("# Adding missing information")
modified_counts = subset_df.copy(deep=True)
for col in col_list:
    modified_counts.loc[modified_counts.sample(frac=0.35).index, col] = pd.np.nan

## randomize 
print ("# Shuffling information")
random_counts = modified_counts.apply(np.random.permutation, axis=1, result_type='broadcast')
random_counts[np.isnan(random_counts)] = 0
random_counts['total'] = random_counts.sum(axis=1)

## get frequencies
print ("# Get frequence")
random_freqs = get_freq(random_counts, col_list)

if (args.debug):
    print ('##########')
    print ('Random original Counts')
    print (subset_df)
    print ('##########')
    print ('')

    print ('##########')
    print ('Random original Frequence')
    subset_df['total'] = subset_df.sum(axis=1)
    original_freq = get_freq(subset_df, col_list)
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
print ("# Adjusting to 100 counts")
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
print (subset_df)
print ('##########')
print ('')

## print randomize counts & frequencies
print ('##########')
print ('Counts adjusted')
print (new_random)
print ('##########')
print ('')

print ("Printing frequencies in table: " + args.out)
#print (df_miRNA)
new_random.to_csv(args.out + ".csv", ',')
