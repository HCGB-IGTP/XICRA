#!/usr/bin/python
import csv
import sys
from collections import defaultdict

"""
This code is copy from PILFER software. See copyright and License details in https://github.com/rishavray/PILFER
Original code: June 2018
https://github.com/rishavray/PILFER/blob/master/tools/annotate_sam.py

Modifications: November 2021
- Add some comments to understand code and clarify it.
- Add to read fasta file as dict and not as list as provided. 

This particular scripts uses input in sam format and a file with gold piRNA sequences to identify putative and well-knonw piRNAs.

"""

def ReverseComplement(seq):
    seq_dict = {'A':'T','T':'A','G':'C','C':'G','N':'N'}
    return "".join([seq_dict[base] for base in reversed(seq)])

### Original code
#seq_id = []
#with open(sys.argv[1],"rb") as listin:
#    for seq in listin:
#        seq_id.append(seq.strip())
#seq_id = list(set(seq_id))

## Modified
from XICRA.scripts import BAMtoPILFER

## read fasta file and save results in list
seq_dict = BAMtoPILFER.get_known_piRNA(sys.argv[1], Debug=False)
seq_id = list(seq_dict.keys())

## parse STDIN sam
mini = 26
maxi = 33

csvin = csv.reader(sys.stdin, delimiter="\t")
csvout = csv.writer(sys.stdout, delimiter="\t")
for row in csvin:
    if not row[0][0] == "@":
        f = row[0].split(":")
        row[0] = f[1]
        seq = row[9]
        if (int(row[1]) & (0x10)):
            seq = ReverseComplement(seq)
        
        ## Check
        if seq in seq_id:
            row.append("XP:Z:PI")
        elif len(row[9])>=mini and len(row[9])<=maxi and int(f[0])>=100:
            row.append("XP:Z:PU")

        row.append("XC:i:"+f[0])
    csvout.writerow(row)
