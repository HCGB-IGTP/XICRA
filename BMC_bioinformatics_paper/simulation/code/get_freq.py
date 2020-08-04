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

## import my modules
from HCGB import functions
from XICRA.scripts import reads2tabular

## dictionary
#freq_fasta = defaultdict(int)

## ARGV
if len (sys.argv) < 4:
    print ("\nUsage:")
    print ("python %s file format out_file\n" %os.path.realpath(__file__))
    exit()

### get line count
n=0
if (sys.argv[2] == 'fastq'):
    n=4
elif (sys.argv[2] == 'fasta'):
    n=2
else:
    print ("** ERROR: please provide a valid format tag: fastq/fasta")
    exit()

# ----------------------------------
# Simulated isomiR categories: ##
# ----------------------------------
# FA <- Five add
# FS <- Five del
# NT <- Non-template
# SR <- SNP Rest
# SS <- SNP Seed
# TA <- Three add
# TS <- Three del
# CN <- Canonical/mature

## dictionary
freq_variants = defaultdict(int)
freq_variants_miRNA = {}

## read file    
with open(sys.argv[1], 'r') as fh:
    lines = []
    for line in fh:
        lines.append(line.rstrip())
        if len(lines) == n:
            record = reads2tabular.process_fasta(lines)
            #print ("##")
            #sys.stderr.write("Record: %s\n" % (str(record)))
            
            # re-init
            lines = []
            
            ## discard: 
            # e.g. >hsa-mir-518f-5p::>hsa-mir-520c-5p|>hsa-mir-526a-5p|>hsa-mir-518d-5p|TS-7511
            # e.g. >hsa-mir-548h-3p::>hsa-mir-548z|TS-5966
            if (re.search('.*::>.*', record['name'])):
                continue

            ## parse the others
            list_split = record['name'].split('::')

            #print (list_split)
            if (sys.argv[2]=='fastq'):
                miRNA = list_split[0].replace('@', '') 
            elif (sys.argv[2]=='fasta'):
                miRNA = list_split[0].replace('>', '')

            variant_list = list_split[1].split('-')
            variant_type = variant_list[0]
            
            ## General
            freq_variants[variant_type] += 1
            
            ## init for miRNA
            if miRNA not in freq_variants_miRNA:
                freq_variants_miRNA[miRNA] = defaultdict(int)
                
            ## add variant type
            freq_variants_miRNA[miRNA][variant_type] += 1
            
##
#for k in sorted (freq_variants.keys()):
#    print("%s\t%s" % (k, freq_variants[k]))

## create pandas
df_miRNA = pd.DataFrame(columns=('FA','FS', 'NT', 'SR', 'SS', 'TA', 'TS', 'CN'))
for miRNA in sorted (freq_variants_miRNA.keys()):
    for k in sorted (freq_variants_miRNA[miRNA].keys()):
        df_miRNA.loc[miRNA, k] = int(freq_variants_miRNA[miRNA][k])

#print (df_miRNA)
df_miRNA.to_csv(sys.argv[3], ',')
    
