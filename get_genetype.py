#usr/bin/env python
import time
import io
import os
import re
import sys
from io import open
from sys import argv
import subprocess

## useful imports

gfffile = argv[1]

## get important lines
parse_gff = "tmp.gff3"
cmd = "awk '{ if ($3 ==\"gene\") {print $0} }' %s > %s" %(gfffile, parse_gff)
#print (cmd)

try:
	subprocess.check_output(cmd, shell = True)
except subprocess.CalledProcessError as err:
	print ('')

## parse
parse_file = open(parse_gff)
text = parse_file.read()
lines = text.splitlines()
	
for line in lines:
	if not line.startswith('#'):
		#print ('## ',line)
		gene = line.split('\t')[-1].split(';')[0].split("ID=")[1]
		gene_type = line.split('\t')[-1].split(';')[2].split("gene_type=")[1]
		print ("%s\t%s" %(gene, gene_type))
