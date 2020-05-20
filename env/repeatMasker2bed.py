#usr/bin/env python
import time
import io
import os
import re
import sys
from io import open
from sys import argv
import pandas as pd

## ARGV
if len (sys.argv) < 4:
	print ("\nUsage:")
	print ("python3 %s repeatMasker info_sequence_names folder\n" %os.path.abspath(argv[0]))
	exit()

repeatMasker_file = argv[1]
conversion_name = argv[2]
folder = argv[3]

df_SeqName = pd.read_csv(conversion_name, sep="\t", header=None, index_col=0, squeeze=True).to_dict()
print ("+ Parsing information provided for sequence names:")

###
split_name = os.path.splitext( os.path.basename(repeatMasker_file) )
repeatmasker_bed = folder + '/' + split_name[0] + '.bed'
# Open file OUT
output_file = open(repeatmasker_bed, 'w')
# Open file IN
fileHandler = open (repeatMasker_file, "r")
while True:
	# Get next line from file
	line = fileHandler.readline()
   	# If line is empty then end of file reached
	if not line :
		break;
	line = line.strip()
	if not line:
		continue
	if line.startswith('SW'):
		continue
	elif line.startswith('score'):
		continue
	else:
		line = re.sub('\s+', '\t', line) ## replace spaces for tabs
		seqID = line.split('\t')[4]
		typeRepeat = line.split('\t')[10]
		if seqID in df_SeqName:
			if typeRepeat != 'Simple_repeat':
				chR = df_SeqName[seqID]	
				start = line.split('\t')[5]
				end = line.split('\t')[6]
				repeatID = line.split('\t')[9]
				score = line.split('\t')[0]
	
				## strand
				strand = '+'
				if line.split('\t')[7] == 'C':
					strand = '-'
				
				string2write = '%s\t%s\t%s\t%s==%s\t%s\t%s\n' %(chR, start, end, repeatID, typeRepeat, score, strand )
				output_file.write(string2write)

## close and finish		
output_file.close()
fileHandler.close()

