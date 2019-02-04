#usr/bin/env python

## useful imports
import time
import io
import os
import re
import sys
from sys import argv
import subprocess

bam_file = os.path.abspath(argv[1])
bedtools_exe = argv[2]
bed_file = os.path.basename(bam_file)

## generate bed file with bedtools bamtobed -i bam_file
cmd_bedtools = "$s bamtobed -i %s > bed_file" %(bedtools_exe, bam_file, bed_file)
try:
	data = cmd_bedtools.result()
except Exception as exc:
	print ('***ERROR:')
	print (cmd_bedtools)
	print('bedtools command generated an exception: %s' %exc)

## paste Aligned.sortedByCoord.out.bed Aligned.sortedByCoord.out.sam | awk -v "OFS=\t" '{print $1, $2, $3, $16, $6}' 


#bedinfo = open(bed_file)
#bedinfo_text = bedinfo.read()
#bedinfo_lines = bedinfo_text.splitlines()
#for line in bedinfo_lines:
#	if not line.startswith('#'):
#		ID = line.split('\t')[3] 
#		
#		if ID in record:
#			print (record.description)
