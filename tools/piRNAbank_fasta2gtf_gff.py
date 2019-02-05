#usr/bin/env python

## useful imports
import time
import io
import os
import re
import sys
from sys import argv
from Bio import SeqIO
from Bio import Seq

input_file = os.path.abspath(argv[1])
name = argv[2]

#print (input_file)
file_gff3 = name + '.gff3'
file_gtf = name + '.gtf'
output_GFF = open(file_gff3, 'w')
output_GTF = open(file_gtf, 'w')

output_GFF.write("##gff-version 3\n")
output_GFF.write("#description: Homo sapiens piRNAbank data formatted from fasta file\n")
output_GFF.write("#provider: piRNAbank\n")
output_GFF.write("#contact: jsanchez\n")
output_GFF.write("#format: gff3\n")
output_GFF.write("#date: 2019-02-01\n")

output_GTF.write("#description: Homo sapiens piRNAbank data formatted from fasta file\n")
output_GTF.write("#provider: piRNAbank\n")
output_GTF.write("#contact: jsanchez\n")
output_GTF.write("##format: gtf\n")
output_GTF.write("#date: 2019-02-01\n")

file_DNA_fasta = name + '.fasta'
output_FASTA = open(file_DNA_fasta, 'w')

for record in SeqIO.parse(input_file, "fasta"):
	#print ("--")
	#print (record.description)

	all_id = record.description
	piRNA_product = all_id.split('|')[0]

	piRNA_ID = ()
	
	gb_search = re.search(".*\|gb\|.*", all_id)
	if gb_search:
		piRNA_ID = all_id.split('|')[2]
	else:
		piRNA_ID = "NULL"
	
	#print (record.seq)
	## transcribe back
	dna_seq = Seq.back_transcribe(record.seq)

	head = ">" + all_id + "\n"
	output_FASTA.write(head)
	output_FASTA.write(str(dna_seq))
	output_FASTA.write("\n")
	
	## parse and generate 
	Chr = all_id.split(':')[1]
	Start = all_id.split(':')[2]
	End = all_id.split(':')[3]
	Strand = all_id.split(':')[-1]

	## GTF example	
	code_strand = ()
	if (Strand == 'Plus'):
		code_strand = '+'
	elif (Strand == 'Minus'):
		code_strand = '-'
	
	
	##	UniSp2  exiqon  spike   1       22      .       +       .       gene_id "UniSp2_rna"; transcript_id "UniSp2_rna"; gene_type "UniSpike"; gene_name "NULL"; transcript_type
	common = "chr" + Chr + "\tpiRNABank\tgene\t" + Start + "\t" + End + "\t.\t" + code_strand + "\t.\t"
	common2 = "chr" + Chr + "\tpiRNABank\ttranscript\t" + Start + "\t" + End + "\t.\t" + code_strand + "\t.\t"
	common3 = "chr" + Chr + "\tpiRNABank\texon\t" + Start + "\t" + End + "\t.\t" + code_strand + "\t.\t"
	
	gff3 = "ID=%s;gene_id=%s;transcript_id=%s;gene_type=piRNA;gene_name=%s;transcript_type=%s;transcript_type=piRNA;level=NULL\n" %(piRNA_product, piRNA_product, piRNA_product, piRNA_ID, piRNA_product)
	gtf ="gene_id \"%s\"; transcript_id \"%s\"; gene_type \"piRNA\"; gene_name \"%s\"; transcript_type \"%s\"; transcript_type \"piRNA\"; level \"NULL\"\n" %(piRNA_product, piRNA_product, piRNA_ID, piRNA_product)

	GFF3 = common + gff3 + common2 + gff3 + common3 + gff3
	GTF = common + gtf + common2 + gtf + common3 + gtf
	
	output_GFF.write(GFF3)
	output_GTF.write(GTF)
	
		
output_GFF.close()
output_GTF.close()
output_FASTA.close()
		
		
		
	
	
	
