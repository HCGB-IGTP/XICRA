#!/usr/bin/env python3
############################################################
## Jose F. Sanchez                                        ##
## Copyright (C) 2019-2020 Lauro Sumoy Lab, IGTP, Spain   ##
############################################################
"""
Get isomiRs from fasta file given a freqs table
"""
## import useful modules
import os
import sys
import re
import time
from io import open
import pandas as pd
from collections import defaultdict
import argparse
import numpy as np

## import my modules
from HCGB import functions
from XICRA.scripts import reads2tabular

#####################################################
parser = argparse.ArgumentParser(prog='get_isomiRs.py', formatter_class=argparse.RawDescriptionHelpFormatter, 
                                 description='''

get_isomiRs.py: Given a frequencies table and fasta file select isomiRs

Version: 0.1
License: GPLv3

USAGE: python get_isomiRs.py --freq table.freq.csv --out out_name --fasta fasta_file [--debug] 
''', epilog="Original code: JFSanchezHerrero")

#####################################################
parser.add_argument('-f', '--freq', action='store', help='Table with original variant frequencies to modify', required=True)
parser.add_argument('-o', '--out', action='store', help='Output names', required=True)
parser.add_argument('--fasta', action='store', help='Fasta file', required=True)
parser.add_argument('--debug', action='store_true', default=False, help='Developer messages')
args = parser.parse_args()
#####################################################

print ("# Read frequency table")
frequencies_miRNA = functions.main_functions.get_data(args.freq, ',', 'index_col=0')
col_list = list(frequencies_miRNA) ## get columns

print ("# Selecting variants from file:" + args.fasta)
print ("# Printing isomiRs sequences in fasta: " + args.out + '.fasta')

## new df
isomiRs_seqs = pd.DataFrame(0, index=frequencies_miRNA.index, columns=col_list)

## read file    
with open(args.out + '.fasta', 'w') as outfh:
    with open(args.fasta, 'r') as fh:
        lines = []
        for line in fh:
            lines.append(line.rstrip())
            if len(lines) == 2:
                record = reads2tabular.process_fasta(lines)
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
                miRNA = list_split[0].replace('>', '') 
                variant_list = list_split[1].split('-')
                variant_type = variant_list[0]
                    
                if miRNA in frequencies_miRNA.index:
                    count_isomiRs = int(frequencies_miRNA.loc[miRNA, variant_type])
                    if (count_isomiRs > 0):
                        
                        if (args.debug):
                            sys.stderr.write("Record: %s\n" % (str(record)))
                            print ("miRNA: " + miRNA)
                            print ("variant: " + variant_type)
                            print ("count: " + str(count_isomiRs))
                            
                        ## Variant already satisfied
                        if (isomiRs_seqs.loc[miRNA, variant_type] != 0):
                            continue
                        
                        ## right now it saves the last 
                        isomiRs_seqs.loc[miRNA, variant_type] = record['name'].replace('>', '') + 'x' + str(count_isomiRs)
    
                        ## print a singlee isomiR for each variant as many times as desired
                        for i in range(0, count_isomiRs):
                            ## string to write: miRNA x times
                            string = record['name'] + '_' + str(i) + '\n' + record['sequence'] + "\n"
                            outfh.write(string)
                    else:
                        isomiRs_seqs.loc[miRNA, variant_type] = ""

##
outfh.close()

print ("# Save frequencies in table: " + args.out + '.csv')
isomiRs_seqs.to_csv(args.out + ".csv", ',')
                    