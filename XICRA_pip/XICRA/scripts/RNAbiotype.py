#!/usr/bin/env python3
############################################################
## Jose F. Sanchez                                        ##
## Copyright (C) 2019-2020 Lauro Sumoy Lab, IGTP, Spain   ##
############################################################
## useful imports
import time
import io
import os
import re
import sys
from sys import argv
import subprocess

## import my modules
from XICRA.config import set_config

## import HCGB
from HCGB.functions import system_call_functions, main_functions
from HCGB.functions import files_functions
##from HCGB.functions import math_functions

## plots
import pandas as pd
import matplotlib

matplotlib.use('agg')
import matplotlib.pyplot as plt
from pandas.plotting import table

#####################
def help_info():
	'''Provide information on RNA biotype analysis'''
	print ("** TODO: Add description on RNAbiotype analysis and details")
	exit()

#####################
def percentage(percent, whole):
	if (percent == "0.00%"):
		return 0

	percent_search = re.search(r"(.*)%", percent)
	if percent_search:
		percent_int = float(percent_search.group(1))
	value = (percent_int * whole) / 100.0
	return int(value)

#####################
def featureCount_call(featureCount_exe, path, gtf_file, bam_file, name, threads, Debug):
	
	## folder for results
	if not os.path.isdir(path):
		files_functions.create_folder(path)

	out_file = os.path.join(path, 'featureCount.out')
	logfile = os.path.join(path, name + '_RNAbiotype.log')

	## debugging messages
	if Debug:
		print ("** DEBUG:")
		print ("featureCounts system call for sample: " + name)
		print ("out_file: " + out_file)
		print ("logfile: " + logfile)

	## send command for feature count
	cmd_featureCount = ('%s --largestOverlap -T %s -p -t exon -g transcript_biotype -a %s -o %s %s 2> %s' %(
		featureCount_exe, threads, gtf_file, out_file, bam_file, logfile)
	)
	## do not count multi mapping, assign only to one feature and in case of ambiguity assign to longest overlap
	#Before: cmd_featureCount = ('%s -M -O -T %s -p -t exon -g transcript_biotype -a %s -o %s %s 2> %s' %(featureCount_exe, threads, gtf_file, out_file, bam_file, logfile))
		
	## system call
	cmd_featureCount_code = system_call_functions.system_call(cmd_featureCount, False, True)
	if not cmd_featureCount_code:
		print("** ERROR: featureCount failed for sample " + name)
		exit()
	
	(extended_Stats, RNAbiotypes_stats) = parse_featureCount(out_file, path, name, bam_file, Debug)

	## debugging messages
	if Debug:
		print ("** DEBUG:")
		print ("extended_Stats: " + extended_Stats)
		print ("RNAbiotypes_stats: " + RNAbiotypes_stats)
		
	return (extended_Stats, RNAbiotypes_stats)

#######################################################################
def parse_featureCount(out_file, path, name, bam_file, Debug):
	"""
	Parses featureCount results for RNAbiotype analysis.
	
	:param out_file: Name provided to featureCount for output results.
	:param path:
	:param name:
	
	
	"""
	
	## debugging messages
	if Debug:
		print ("** DEBUG:")
		print ("Parse results for sample: " + name)
		
		
	## parse results
	out_tsv_file_name = out_file + '.tsv'
	out_tsv_file = open(out_tsv_file_name, 'w')
	RNA_biotypes_file_name = os.path.join(path, name + '_RNAbiotype.tsv')
	RNA_biotypes_file = open(RNA_biotypes_file_name, 'w')
	tRNA_count = 0
	
	##########################################
	### read count file
	##########################################
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
	RNA_biotypes_file.close()
			
	##########################################
	### read summary count file
	##########################################
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

	##########################################
	## get mapping statistics according to mapping software
	##########################################
	count_multi = 0
	count_unmap = 0
	mapping_folder = os.path.dirname(bam_file)
	mapping_stats = mapping_folder + '/Log.final.out'
	
	## -------------------------------- ##
	### STAR mapping		
	## -------------------------------- ##
	if files_functions.is_non_zero_file(mapping_stats):
		## debugging messages
		if Debug:
			print ("** DEBUG:")
			print ("STAR mapping available for sample: " + name)
			print ("mapping_folder: " + mapping_folder)

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
				#count_reads = math_functions.percentage(perc_tmp, total_input_reads)
				count_unmap = count_unmap + count_reads
	else:

		## -------------------------------- ##
		## tophat
		## -------------------------------- ##

		mapping_stats = mapping_folder + '/align_summary.txt' 
		count_map = 0
		total_input_reads = 0
		
		if files_functions.is_non_zero_file(mapping_stats):
			## debugging messages
			if Debug:
				print ("** DEBUG:")
				print ("tophat mapping available for sample: " + name)
				print ("mapping_folder: " + mapping_folder)
			
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
	string2write_unmap = "unmapped\t%s\n" %count_unmap
	out_tsv_file.write(string2write_unmap)
	
	## close files
	out_tsv_file.close()

	return(out_tsv_file_name, RNA_biotypes_file_name)

#######################################################################
def RNAbiotype_module_call(samples_dict, output_dict, gtf_file, threads, Debug):
	"""
	Create RNAbiotype analysis for each sample and create summary plots
	
	:param samples_dict: Dictionary containing sample IDs as keys and bam files and folder for results as values
	:param gtf_file: Gene annotation file for the reference genome used.
	:param threads: Number of threads to use.
	:param Debug: True/False for debugging messages
	"""
	
	## get bin
	featureCount_exe = set_config.get_exe('featureCounts')

	## loop dictionary
	for sample in samples_dict.items():
		print (samples_dict[sample])
		print (output_dict[sample])
		
		#(extended_Stats, RNAbiotypes_stats) = featureCount_call(featureCount_exe, folder, gtf_file, bam_file, name, threads, Debug)
		#pie_plot_results(RNAbiotypes_stats, Debug)


	return()

############################################################
def generate_matrix(dict_files):
	######################################################
	### For a dictionary containing names as keys and files as values, 
	###	generates a count matrix with the classification from featurecounts
	###		
	### It returns a dataframe containing for each index generated
	### count values for each sample (sample_name) in columns
	###
	#########################################################
	all_data = pd.DataFrame()
	for key,values in dict_files.items():
		print ('+ Reading information from sample: ', key)	
		data = pd.read_csv(values, sep='\t', header=None, names=['RNAbiotypes', key])
		
		## skip if file is empty
		if data.empty:
			continue

		## get info, generate unique name and merge for samples
		data = data.set_index('RNAbiotypes')
		all_data = pd.concat([all_data, data], axis=1, sort=True)
	
	return (all_data)
	##

#######################################################################
def pie_plot_results(RNAbiotypes_stats_file, Debug):
	# PLOT and SHOW results
	RNAbiotypes_stats = main_functions.get_data(RNAbiotypes_stats_file, '\t', 'index_col=0')

	# create plot
	plt.figure(figsize=(16,8))
	df_genetype_2 = pd.DataFrame({'Type':RNAbiotypes_stats[0], 'Read_Count':RNAbiotypes_stats[1]}).sort_values(by=['Read_Count'])

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


#######################################################################
def main():
	
	## ARGV
	if len (sys.argv) < 6:
		print ("\nUsage:")
		print ("python3 %s bam_file folder gtf_file threads name featureCount_bin\n" %os.path.realpath(__file__))
		exit()
	
	bam_file = os.path.abspath(argv[1])
	folder = os.path.abspath(argv[2])
	gtf_file = os.path.abspath(argv[3])
	threads = argv[4]
	name = argv[5]
	featureCount_exe = argv[6]

	## Debug
	Debug=True
	
	## variables
	(extended_Stats_file, RNAbiotypes_stats_file) = featureCount_call(featureCount_exe, folder, gtf_file, bam_file, name, threads, Debug)
	pie_plot_results(RNAbiotypes_stats_file, Debug)
	
	
######
if __name__== "__main__":
	main()
	
