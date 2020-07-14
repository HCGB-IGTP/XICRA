#!/usr/bin/env python3

## useful imports
import time
import io
import os
import re
import sys
from sys import argv
from collections import defaultdict

## shell: grep -v '>' file.fasta | sort | uniq -c | awk '{print $2"\t"$1}' > tab_info.txt
## source: https://www.biostars.org/p/317524/

### 
def process_fasta(lines):
    ks = ['name', 'sequence']
    return {k: v for k, v in zip(ks, lines)}

### 
def process_fastq(lines):
    ks = ['name', 'sequence', 'optional', 'quality']
    return {k: v for k, v in zip(ks, lines)}

#####
def reads2tabular(fastq_file, out):
    
    ## dictionary
    freq_fasta = defaultdict(int)
    
    ## read fastq    
    n = 4
    with open(fastq_file, 'r') as fh:
        lines = []
        for line in fh:
            lines.append(line.rstrip())
            if len(lines) == n:
                record = process_fastq(lines)
                #sys.stderr.write("Record: %s\n" % (str(record)))
                lines = []
                
                ## add sequences & count
                freq_fasta[record['sequence']] += 1

    ## print in file
    with open(out,'w') as file:
        for k in sorted (freq_fasta.keys()):
            file.write("%s\t%s\n" % (k, freq_fasta[k]))
    
    return(freq_fasta)
            
######
def help_options():
    print ("\nUSAGE: python %s fastq out_file\n"  %os.path.realpath(__file__))

######
def main():
    ## this code runs when call as a single script

    ## control if options provided or help
    if len(sys.argv) > 1:
        print ("")
    else:
        help_options()
        exit()        
    
    reads2tabular(argv[1], argv[2])
        
######
if __name__== "__main__":
    main()
    
    
