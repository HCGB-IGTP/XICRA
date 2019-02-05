#usr/bin/env python

## useful imports
import time
import io
import os
import re
import sys
from sys import argv
import subprocess

## ARGV
if len (sys.argv) < 4:
	print ("\nUsage:")
	print ("python3 %s bam_file bedtools_bin samtools_bin logfile\n" %argv[0])
	exit()

bam_file = os.path.abspath(argv[1])
bedtools_exe = argv[2]
samtools_exe = argv[3]
logFile = argv[4]

# start
output_file = open(logFile, 'a')
output_file.write("\nConvert BAM to Pilfer Input file:\n")	

## Variables
dirname_name = os.path.dirname(bam_file)
split_name = os.path.splitext( os.path.basename(bam_file) )
bed_file = dirname_name + '/' + split_name[0] + '.bed'
sam_file = dirname_name + '/' + split_name[0] +  '.sam'
pilfer_tmp = dirname_name + '/' + split_name[0] + '.tmp.pilfer.bed'
pilfer_file = dirname_name + '/' + split_name[0] + '.pilfer.bed'

## START
print ("\n+ Converting BAM file into PILFER input file")

## generate bed file with bedtools bamtobed -i bam_file
if (os.path.isfile(bed_file)):
	print ("\t+ File %s already exists" %bed_file)
else:
	cmd_bedtools = "%s bamtobed -i %s > %s" %(bedtools_exe, bam_file, bed_file)
	output_file.write(cmd_bedtools)
	output_file.write("\n")		
	try:
		subprocess.check_output(cmd_bedtools, shell = True)
	except Exception as exc:
		print ('***ERROR:')
		print (cmd_bedtools)
		print('bedtools command generated an exception: %s' %exc)
		exit()
		
## generate samtools
if (os.path.isfile(sam_file)):
	print ("\t+ File %s already exists" %sam_file)
else:
	cmd_samtools = "%s view %s > %s" %(samtools_exe, bam_file, sam_file)
	output_file.write(cmd_samtools)
	output_file.write("\n")		
	try:
		subprocess.check_output(cmd_samtools, shell = True)
	except Exception as exc:
		print ('***ERROR:')
		print (cmd_samtools)
		print('samtools view command generated an exception: %s' %exc)
		exit()

## generate paste filter tmp file
if (os.path.isfile(pilfer_tmp)):
	print ("\t+ File %s already exists" %pilfer_tmp)
else:
	## paste Aligned.sortedByCoord.out.bed Aligned.sortedByCoord.out.sam | awk -v "OFS=\t" '{print $1, $2, $3, $16, $6}' 
	cmd_paste = "paste %s %s | awk -v \"OFS=\t\" \'{print $1, $2, $3, $16, $6}\' > %s" %(bed_file, sam_file, pilfer_tmp)
	output_file.write(cmd_paste)
	output_file.write("\n")		
	try:
		subprocess.check_output(cmd_paste, shell = True)
	except Exception as exc:
		print ('***ERROR:')
		print (cmd_paste)
		print('paste bed sam command generated an exception: %s' %exc)
		exit()

## parse pilfer tmp file
counter = 1
previous_line = ()

# Open file OUT
output_file = open(pilfer_file, 'w')

# Open file IN
fileHandler = open (pilfer_tmp, "r")
while True:
	# Get next line from file
	line = fileHandler.readline().strip()
   	# If line is empty then end of file reached
	if not line :
		break;

	seq = line.split('\t')[3]
	real_seq = seq.split('::PI')
	seq_len = len(str(real_seq[0]))

	## Discard smaller
	if (previous_line):
		if (previous_line == line):
			line = previous_line
			counter += 1
		else:
			line_split = previous_line.split('\t')
			output_file.write('%s\t%s\t%s\t%s::PI\t%s\t%s\n' %(line_split[0], line_split[1], line_split[2], line_split[3], counter, line_split[4]))

	#counter += 1
	while True:
		#get next line
		next_line = fileHandler.readline().strip()
	
		if (next_line == line):
			counter += 1
		else:
			line_split = line.split('\t')
			output_file.write('%s\t%s\t%s\t%s::PI\t%s\t%s\n' %(line_split[0], line_split[1], line_split[2], line_split[3], counter, line_split[4]))

			previous_line = next_line
			counter = 1
			break;
 
# Close Close    
fileHandler.close()
output_file.close()
