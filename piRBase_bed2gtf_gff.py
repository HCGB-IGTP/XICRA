#usr/bin/env python

## useful imports
import time
import io
import os
import re
import sys
from sys import argv

input_file = os.path.abspath(argv[1])
name = argv[2]

#print (input_file)
file_gff3 = name + '.gff3'
file_gtf = name + '.gtf'
output_GFF = open(file_gff3, 'w')
output_GTF = open(file_gtf, 'w')

output_GFF.write("##gff-version 3\n")
output_GFF.write("#description: Homo sapiens piRBase data formatted from BED file\n")
output_GFF.write("#provider: piRBase\n")
output_GFF.write("#contact: jsanchez\n")
output_GFF.write("#format: gff3\n")
output_GFF.write("#date: 2019-02-04\n")

output_GTF.write("#description: Homo sapiens piRBase data formatted from fasta file\n")
output_GTF.write("#provider: piRBase\n")
output_GTF.write("#contact: jsanchez\n")
output_GTF.write("##format: gtf\n")
output_GTF.write("#date: 2019-02-04\n")

bedinfo = open(input_file)
bedinfo_text = bedinfo.read()
bedinfo_lines = bedinfo_text.splitlines()
for line in bedinfo_lines:
	if not line.startswith('#'):

		## parse and generate 
		Chr = line.split('\t')[0]
		Start = line.split('\t')[1]
		End = line.split('\t')[2]
		strand = line.split('\t')[5]
		ID= line.split('\t')[3]

		common = Chr + "\tpiRBase\tgene\t" + Start + "\t" + End + "\t.\t" + strand + "\t.\t"
		common2 = Chr + "\tpiRBase\ttranscript\t" + Start + "\t" + End + "\t.\t" + strand + "\t.\t"
		common3 = Chr + "\tpiRBase\texon\t" + Start + "\t" + End + "\t.\t" + strand + "\t.\t"
	
		gff3 = "ID=%s;gene_id=%s;transcript_id=%s;gene_type=piRNA;gene_name=%s;transcript_type=%s;transcript_type=piRNA;level=NULL\n" %(ID, ID, ID, ID, ID)
		gtf ="gene_id \"%s\"; transcript_id \"%s\"; gene_type \"piRNA\"; gene_name \"%s\"; transcript_type \"%s\"; transcript_type \"piRNA\"; level \"NULL\"\n" %(ID, ID, ID, ID)

		GFF3 = common + gff3 + common2 + gff3 + common3 + gff3
		GTF = common + gtf + common2 + gtf + common3 + gtf
	
		output_GFF.write(GFF3)
		output_GTF.write(GTF)
		
output_GFF.close()
output_GTF.close()
	
	
