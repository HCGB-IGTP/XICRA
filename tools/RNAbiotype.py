#usr/bin/env python

## useful imports
import time
import io
import os
import re
import sys
from sys import argv
import subprocess

## plots
import pandas as pd
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from pandas.plotting import table

## ARGV
if len (sys.argv) < 6:
	print ("\nUsage:")
	print ("python3 %s bam_file folder gtf_file featureCount_bin logFile threads\n" %os.path.realpath(__file__))
	exit()

bam_file = os.path.abspath(argv[1])
folder = argv[2]
gtf_file = argv[3]
featureCount_exe = argv[4]
logFile = argv[5]
threads = argv[6]

## start
output_file = open(logFile, 'a')

## variables
#featureCount_exe = config['EXECUTABLES']['featureCount_exe']
logfile = folder + '/featureCount_logfile.txt'
out_file = folder + '/featureCount.out'
out_tsv = folder + '/featureCount.out.tsv'
RNA_biotypes = folder + '/RNAbiotypes.tsv'
name_figure = folder + '/RNAbiotypes.pdf'

if not (os.path.isfile(bam_file)):
	print ('***ERROR: bam file does not exist')
	exit()

else:

	if not (os.path.isfile(out_file)):
		output_file.write("\nParse RNA Biotype results:\n")	
		
		## send command for feature count
		cmd_featureCount = ('%s -M -O -T %s -p -t exon -g transcript_biotype -a %s -o %s %s 2> %s' %(featureCount_exe, threads, gtf_file, out_file, bam_file, logfile))
		
		print (cmd_featureCount)
		output_file.write(cmd_featureCount)
		output_file.write("\n")		
		# send command	
		try:
			subprocess.check_output(str(cmd_featureCount), shell = True)
		except subprocess.CalledProcessError as err:
			print ('***ERROR:')
			print ('[CMD: %s ]' %cmd_featureCount)
			print (err.output)
	else:
		print ("Feature counts are already analyzed...")
		
	## parse count
	## prepare file for plot
	if not (os.path.isfile(name_figure)):
		count_file = open(out_file)
		count_file_text = count_file.read()
		count_file_lines = count_file_text.splitlines()	
		out_tsv_file = open(out_tsv, 'w')

		RNA_biotypes_file = open(RNA_biotypes, 'w')
		tRNA_count = 0
	
		for line in count_file_lines:
			if line.startswith('#'):
				continue
			elif line.startswith('Geneid'):
				continue
			else:
				
				ID = line.split('\t')[0]
				count = int(line.split('\t')[-1])
				string2write_raw = "%s\t%s\n" %(ID, count)
				out_tsv_file.write(string2write_raw)

				tRNA_search = re.search(r".*tRNA", ID)
				if tRNA_search:
					tRNA_count = int(tRNA_count) + int(count)				
				elif (count > 0):
					RNA_biotypes_file.write(string2write_raw)
		
		string2write = "tRNA\t%s\n" %tRNA_count
		RNA_biotypes_file.write(string2write)
		RNA_biotypes_file.close()
		out_tsv_file.close()

		output_file.close()
	
		# PLOT and SHOW results
		## parse results	
		df_genetype = pd.read_csv(RNA_biotypes, sep="\t", header=None)	

		# create plot
		plt.figure(figsize=(16,8))
		df_genetype_2 = pd.DataFrame({'Type':df_genetype[0], 'Read_Count':df_genetype[1]}).sort_values(by=['Read_Count'])

		## get total count
		df_genetype_ReadCount_sum = df_genetype_2['Read_Count'].sum()

		## filter 1% values
		minimun = df_genetype_ReadCount_sum * 0.01
		df_genetype_filter_greater = df_genetype_2[ df_genetype_2['Read_Count'] >= minimun ]
		df_genetype_filter_smaller = df_genetype_2[ df_genetype_2['Read_Count'] < minimun ]
	
		## merge and generate Other class
		df_genetype_filter_smaller_sum = df_genetype_filter_smaller['Read_Count'].sum() ## total filter smaller
		df_genetype_filter_greater2 = df_genetype_filter_greater.append({'Read_Count':df_genetype_filter_smaller_sum, 'Type':'Other'}, ignore_index=True)
	
	
		## Create Plot
		ax1 = plt.subplot(121, aspect='equal')
		df_genetype_filter_greater2.plot.pie(y = 'Read_Count', ax=ax1, autopct='%1.2f%%', shadow=False, labels=df_genetype_filter_greater2['Type'], legend = False)

		# plot table
		ax2 = plt.subplot(122)
		plt.axis('off')
		tbl = ax2.table(cellText=df_genetype_2.values, colLabels=df_genetype_2.columns, loc='center', rowLoc='left', cellLoc='center', colWidths=[0.15, 0.5])
		tbl.auto_set_font_size(True)
		tbl.scale(1.1,1.1)
	
		## generate image
		plt.savefig(name_figure)
		
	else:
		print ("PDF plot is already generated...")
