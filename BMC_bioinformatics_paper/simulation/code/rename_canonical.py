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

## import my modules
from HCGB.functions import fasta_functions

## ARGV
if len (sys.argv) < 2:
    print ("\nUsage:")
    print ("python %s fasta_file out_file\n" %os.path.realpath(__file__))
    exit()
    
## read fastq    
n = 2
with open(sys.argv[1], 'r') as fh:
    with open(sys.argv[2],'w') as file:
            
        lines = []
        for line in fh:
            lines.append(line.rstrip())
            if len(lines) == n:
                record = fasta_functions.process_fasta(lines)
                #print ("##")
                #sys.stderr.write("Record: %s\n" % (str(record)))
                
                lines = []
                ## add tag
                new_name = record['name'].split(' ')[0] + '::CN'
                new_name = new_name.replace('miR', 'mir')
                
                ## convert U->T
                new_seq = record['sequence'].replace('U', 'T')
                
                ## print in file
                file.write("%s\n%s\n" % (new_name, new_seq))

