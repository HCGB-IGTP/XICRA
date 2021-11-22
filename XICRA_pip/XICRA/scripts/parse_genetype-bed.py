#usr/bin/env python
import time
import io
import os
import re
import sys
from io import open
from sys import argv
import subprocess

## ARGV
if len (sys.argv) < 2:
    print ("\nUsage:")
    print ("python3 %s gtffile folder\n" %os.path.realpath(__file__))
    exit()

gtfFile = os.path.abspath(argv[1])
folder = argv[2]

## Variables
dirname_name = os.path.dirname(gtfFile)
split_name = os.path.splitext( os.path.basename(gtfFile) )

# Open OUT
exon_gtf = folder + '/' + split_name[0] + '_exon.gtf' 
output_exon = open(exon_gtf, 'w')
miRNA_gtf = folder + '/' + split_name[0] + '_miRNA.gtf' 
output_miRNA = open(miRNA_gtf, 'w')
ncRNA_bed = folder + '/' + split_name[0] + '_ncRNA.bed' 
output_ncRNA = open(ncRNA_bed, 'w')

# Open file IN
fileHandler = open(gtfFile, "r")
while True:
    # Get next line from file
    line = fileHandler.readline().strip()
       # If line is empty then end of file reached
    if not line :
        break;
    
    if not line.startswith('#'):        
        typeID = line.split('\t')[2]
        
        if (typeID == 'exon'):
            #print (line)        
            output_exon.write(line)
            output_exon.write('\n')

        description = line.split('\t')[-1]
        description_list = line.split('\t')[-1].split(';')
        gene = []
        gene_type = []

        for items in description_list:
            gene_type_search = re.search('gene_type "(.*)"', items)
            gene_name_search = re.search('gene_name "(.*)"', items)
            if gene_type_search:
                gene_type = gene_type_search.group(1)
            elif gene_name_search:
                gene = gene_name_search.group(1)

        gene_type = str(gene_type)
        gene = str(gene)
        
        chR = line.split('\t')[0]    
        start = line.split('\t')[3]
        end = line.split('\t')[4]
        score = line.split('\t')[5]
        strand = line.split('\t')[6]
        source = line.split('\t')[1]
        extra = line.split('\t')[7]
        
        # write
        string2write = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %(chR, start, end, gene, score, strand, source, typeID, extra, description)
        
        if gene_type in ["misc_RNA", "rRNA", "miRNA", "sRNA", "snRNA", "snoRNA", "scRNA", "macro_lncRNA", "scaRNA", "rRNA", "lincRNA"]:
            output_ncRNA.write(string2write)
            
        if (re.search('(.*)tRNA', gene_type)):
            output_ncRNA.write(string2write)        
        ## we are not writing piRNA as we would like to substract that information

        if (gene_type == 'miRNA'): ## print only miRNA gtf for miRNA analysis
            output_miRNA.write(line)
            output_miRNA.write('\n')
            
## close and finish
output_miRNA.close()
output_exon.close()
output_ncRNA.close()
fileHandler.close()            
