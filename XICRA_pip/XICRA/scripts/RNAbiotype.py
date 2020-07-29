#!/usr/bin/env python3

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

######
def percentage(percent, whole):
	if (percent == "0.00%"):
		return 0

	percent_search = re.search(r"(.*)%", percent)
	if percent_search:
		percent_int = float(percent_search.group(1))
	value = (percent_int * whole) / 100.0
	return int(value)

#####################
if not (os.path.isfile(bam_file)):
	print ('***ERROR: bam file does not exist')
	exit()

else:

	if not (os.path.isfile(out_file)):
		output_file.write("\nParse RNA Biotype results:\n")	
		
		## send command for feature count
		cmd_featureCount = ('%s --largestOverlap -T %s -p -t exon -g transcript_biotype -a %s -o %s %s 2> %s' %(featureCount_exe, threads, gtf_file, out_file, bam_file, logfile))
		
		## do not count multi mapping, assign only to one feature and in case of ambiguity assign to longest overlap
		#Before: cmd_featureCount = ('%s -M -O -T %s -p -t exon -g transcript_biotype -a %s -o %s %s 2> %s' %(featureCount_exe, threads, gtf_file, out_file, bam_file, logfile))
		
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

		out_tsv_file = open(out_tsv, 'w')
		RNA_biotypes_file = open(RNA_biotypes, 'w')
		tRNA_count = 0

		### read count file
		count_file = open(out_file)
		count_file_text = count_file.read()
		count_file_lines = count_file_text.splitlines()	

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
		
		## count and summary tRNA
		string2write = "tRNA\t%s\n" %tRNA_count
		RNA_biotypes_file.write(string2write)
				
		### read summary count file
		summary_count_file = open(out_file + '.summary')
		summary_count_file_text = summary_count_file.read()
		summary_count_file_lines = summary_count_file_text.splitlines()	

		for line in summary_count_file_lines:
			if line.startswith('Status'):
				continue
			elif line.startswith('Assigned'):
				continue
			else:
				## adds Unassigned_Ambiguity
				## adds Unassigned_NoFeatures
				ID = line.split('\t')[0]
				count = int(line.split('\t')[-1])

				## skip empty entries
				if count == 0:
					continue
				string2write_raw = "%s\t%s\n" %(ID, count)
				out_tsv_file.write(string2write_raw)
		
		## get mapping statistics according to mapping software
		count_multi = 0
		count_unmap = 0
		mapping_folder = os.path.dirname(bam_file)
		mapping_stats = mapping_folder + '/Log.final.out'
		
		### STAR mapping		
		if os.path.isfile(mapping_stats):
			mapping_stats_file = open(mapping_stats)
			mapping_stats_file_text = mapping_stats_file.read()
			mapping_stats_file_lines = mapping_stats_file_text.splitlines()	
	
			for line in mapping_stats_file_lines:
				multi_search = re.search(r".*Number of reads mapped to", line)
				unmap_search = re.search(r".*unmapped.*", line)
				input_search = re.search(r".*input reads.*", line)
			
				if input_search:
					total_input_reads = int(line.split('\t')[-1])

				if multi_search:
					count_tmp = int(line.split('\t')[-1])
					count_multi = count_multi + count_tmp

				elif unmap_search:
					perc_tmp = line.split('\t')[-1]
					count_reads = percentage(perc_tmp, total_input_reads)
					count_unmap = count_unmap + count_reads
		else:
			mapping_stats = mapping_folder + '/align_summary.txt' 
			count_map = 0
			total_input_reads = 0
			## tophat
			if os.path.isfile(mapping_stats):
				mapping_stats_file = open(mapping_stats)
				mapping_stats_file_text = mapping_stats_file.read()
				mapping_stats_file_lines = mapping_stats_file_text.splitlines()	
	
				for line in mapping_stats_file_lines:
					map_search2 = re.search(r"Aligned.*\:\s+(\d+).*", line)
					input_search2 = re.search(r".*Input.*\:\s+(\d+).*", line)
					if input_search2:
						total_input_reads = input_search2.group(1)
					if map_search2:
						count_map = map_search2.group(1)
		
				####
				count_unmap = int(total_input_reads) - int(count_map)

			else:
				## other
				print ("Neither tophat or STAR..., no mapping statistics")

		### print mapping stats
		#string2write_multi = "multimapping\t%s\n" %count_multi
		#out_tsv_file.write(string2write_multi)
		
		string2write_unmap = "unmapped\t%s\n" %count_unmap
		out_tsv_file.write(string2write_unmap)
		
		#string2write_total = "total\t%s\n" %total_input_reads
		#out_tsv_file.write(string2write_total)
		
		## close files
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
