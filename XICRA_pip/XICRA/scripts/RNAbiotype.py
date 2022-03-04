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
from termcolor import colored
import concurrent.futures

## import my modules
from XICRA.config import set_config

## import HCGB
from HCGB.functions import system_call_functions, main_functions, time_functions
from HCGB.functions import files_functions, math_functions

## plots
import pandas as pd
import matplotlib
from HCGB.functions.aesthetics_functions import debug_message

matplotlib.use('agg')
import matplotlib.pyplot as plt
from pandas.plotting import table

#####################
def help_info():
	'''Provide information on RNA biotype analysis'''
	print ("** TODO: Add description on RNAbiotype analysis and details")
	exit()

#####################
def biotype_all(featureCount_exe, path, gtf_file, bam_file, name, threads, Debug, allow_multimap, stranded):
	
	## folder for results
	if not os.path.isdir(path):
		files_functions.create_folder(path)

	out_file = os.path.join(path, 'featureCount.out')
	logfile = os.path.join(path, name + '_RNAbiotype.log')

	filename_stamp_all = path + '/.success_all'
	if os.path.isfile(filename_stamp_all):
		stamp = time_functions.read_time_stamp(filename_stamp_all)
		print (colored("\tA previous command generated results on: %s [%s -- %s]" %(stamp, name, 'RNAbiotype'), 'yellow'))
		return()

	else:
		filename_stamp_featureCounts = path + '/.success_featureCounts'
		if os.path.isfile(filename_stamp_featureCounts):
			stamp = time_functions.read_time_stamp(filename_stamp_featureCounts)
			print (colored("\tA previous command generated results on: %s [%s -- %s]" %(stamp, name, 'featureCounts'), 'yellow'))
		else:
            
			## debugging messages
			if Debug:
				print ("** DEBUG:")
				print ("featureCounts system call for sample: " + name)
				print ("out_file: " + out_file)
				print ("logfile: " + logfile)
		
			## send command for feature count
			## Allow multimapping
			if allow_multimap:
				cmd_featureCount = ('%s -s %s -M -O -T %s -p -t exon -g transcript_biotype -a %s -o %s %s 2> %s' %(
					featureCount_exe, stranded, threads, gtf_file, out_file, bam_file, logfile)
				)
			else:
				cmd_featureCount = ('%s -s %s --largestOverlap -T %s -p -t exon -g transcript_biotype -a %s -o %s %s 2> %s' %(
					featureCount_exe, stranded, threads, gtf_file, out_file, bam_file, logfile)
				)
				
				
			## system call
			cmd_featureCount_code = system_call_functions.system_call(cmd_featureCount, False, True)
			if not cmd_featureCount_code:
				print("** ERROR: featureCount failed for sample " + name)
				exit()
				
			## print time stamp
			time_functions.print_time_stamp(filename_stamp_featureCounts)
		
		## parse results
		(extended_Stats_file, RNAbiotypes_stats_file) = parse_featureCount(out_file, path, name, bam_file, Debug)
		
		## debugging messages
		if Debug:
			print ("** DEBUG:")
			print ("extended_Stats: " + extended_Stats_file)
			print (main_functions.get_data(extended_Stats_file, '\t', 'header=None'))
			print ("RNAbiotypes_stats: " + RNAbiotypes_stats_file)
			print (main_functions.get_data(RNAbiotypes_stats_file, '\t', 'header=None'))

	return ()

#######################################################################
def parse_featureCount(out_file, path, name, bam_file, Debug):
	"""
	Parses featureCount results for RNAbiotype analysis.
	
	:param out_file: Name provided to featureCount for output results.
	:param path:
	:param name:
	
	
	"""

	## file names
	out_tsv_file_name = out_file + '.tsv'
	RNA_biotypes_file_name = os.path.join(path, name + '_RNAbiotype.tsv')

	##
	filename_stamp_parse = path + '/.success_parse'
	if os.path.isfile(filename_stamp_parse):
		stamp = time_functions.read_time_stamp(filename_stamp_parse)
		print (colored("\tA previous command generated results on: %s [%s -- %s]" %(stamp, name, 'parse results'), 'yellow'))
	else:
	
		## debugging messages
		if Debug:
			print ("** DEBUG:")
			print ("Parse results for sample: " + name)
			
		## parse results
		out_tsv_file = open(out_tsv_file_name, 'w')
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
		
		## ATTENTION: See example at the end of this file
		
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
				unmap_search = re.search(r".*Number of reads unmapped", line)
				input_search = re.search(r".*input reads.*", line)
			
				if Debug:
					debug_message("mapping_stats_file Line")
					debug_message(line)
			
				if input_search:
					total_input_reads = int(line.split('\t')[-1])
					if Debug:
						print("total_input_reads")
						print(type(total_input_reads))
						print(total_input_reads)
							
				if multi_search:
					count_tmp = int(line.split('\t')[-1])
					count_multi = count_multi + count_tmp
					if Debug:
						print("count_tmp")
						print(type(count_tmp))
						print(count_tmp)
						print("count_multi")
						print(type(count_multi))
						print(count_multi)
	
				if unmap_search:
					count_reads_tmp = int(line.split('\t')[-1])
					#count_reads = math_functions.percentage(perc_tmp, total_input_reads)
					count_unmap = count_unmap + count_reads_tmp
					if Debug:
						print("perc_tmp")
						print(type(count_reads_tmp))
						print(count_reads_tmp)
						
						print("count_unmap")
						print(type(count_unmap))
						print(count_unmap)
		
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
		summary_count_file.close()
		mapping_stats_file.close()
		count_file.close()
		## print timestamp
		time_functions.print_time_stamp(filename_stamp_parse)

	return(out_tsv_file_name, RNA_biotypes_file_name)

#######################################################################
def RNAbiotype_module_call(samples_dict, output_dict, gtf_file, Debug, max_workers_int, threads_job, multimapping, stranded):
	"""
	Create RNAbiotype analysis for each sample and create summary plots
	
	:param samples_dict: Dictionary containing sample IDs as keys and bam files as values
	:param output_dict: Dictionary containing sample IDs as keys and output folder as values
	:param gtf_file: Gene annotation file for the reference genome used.
	:param threads: Number of threads to use.
	:param Debug: True/False for debugging messages
	"""
	
	## get bin
	featureCount_exe = set_config.get_exe('featureCounts')

	## send for each sample
	with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers_int) as executor:
		commandsSent = { executor.submit(biotype_all, featureCount_exe, 
										output_dict[sample], gtf_file, bam_files, 
										sample, threads_job, Debug, multimapping, stranded): sample for sample, bam_files in samples_dict.items() }
	
		for cmd2 in concurrent.futures.as_completed(commandsSent):
			details = commandsSent[cmd2]
			try:
				data = cmd2.result()
			except Exception as exc:
				print ('***ERROR:')
				print (cmd2)
				print('%r generated an exception: %s' % (details, exc))

	##
	## plot results
	for name, folder in output_dict.items():
		RNAbiotypes_stats_file = os.path.join(folder, name + '_RNAbiotype.tsv')
		if files_functions.is_non_zero_file(RNAbiotypes_stats_file):
			pie_plot_results(RNAbiotypes_stats_file, name, folder, Debug)
			
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
def pie_plot_results(RNAbiotypes_stats_file, name, folder, Debug):
	
	##
	filename_stamp_plot = folder + '/.success_plot'
	if os.path.isfile(filename_stamp_plot):
		stamp = time_functions.read_time_stamp(filename_stamp_plot)
		print (colored("\tA previous command generated results on: %s [%s -- %s]" %(stamp, name, 'plot results'), 'yellow'))
	else:
		
		# PLOT and SHOW results
		RNAbiotypes_stats = main_functions.get_data(RNAbiotypes_stats_file, '\t', 'header=None')
	
		# create plot
		plt.figure(figsize=(16,8))
		df_genetype_2 = pd.DataFrame({'Type':RNAbiotypes_stats[0], 
									'Count':RNAbiotypes_stats[1]}).sort_values(by=['Count'])
	
		## get total count
		df_genetype_ReadCount_sum = df_genetype_2['Count'].sum()
	
		## filter 1% values
		minimun = df_genetype_ReadCount_sum * 0.01
		df_genetype_filter_greater = df_genetype_2[ df_genetype_2['Count'] >= minimun ]
		df_genetype_filter_smaller = df_genetype_2[ df_genetype_2['Count'] < minimun ]
	
		## create %values
		df_genetype_2['Percentage'] = (df_genetype_2['Count']/df_genetype_ReadCount_sum*100).round(3)
		
		## merge and generate Other class
		df_genetype_filter_smaller_sum = df_genetype_filter_smaller['Count'].sum() ## total filter smaller
		df_genetype_filter_greater2 = df_genetype_filter_greater.append({
			'Count':df_genetype_filter_smaller_sum, 
			'Type':'Other'}, ignore_index=True)
	
		## Create Pie Plot
		ax1 = plt.subplot(121, aspect='equal')
		df_genetype_filter_greater2.plot.pie(
			y = 'Count', 
			ax=ax1, 
			autopct='%1.2f%%', 
			shadow=False, 
			labels=df_genetype_filter_greater2['Type'], 
			legend = False)
	
		# plot table
		ax2 = plt.subplot(122)
		plt.axis('off')
		tbl = ax2.table(
			cellText=df_genetype_2.values, 
			colLabels=df_genetype_2.columns,
			loc='center', rowLoc='left', cellLoc='center', 
			)
		tbl.auto_set_font_size(True)
		#tbl.set_fontsize(12)
		tbl.scale(1.1,1.1)
	
		## set PDF name
		name_figure = os.path.join(folder, name + '_RNAbiotypes.pdf')
	
		## generate image
		plt.savefig(name_figure)		
		plt.close(name_figure)
		plt.close()
		
		## print time stamps
		time_functions.print_time_stamp(filename_stamp_plot)
		filename_stamp_all = folder + '/.success_all'
		time_functions.print_time_stamp(filename_stamp_all)
		
#######################################################################
def main():
	
	## ARGV
	if len (sys.argv) < 6:
		print ("\nUsage:")
		print ("python3 %s bam_file folder gtf_file threads name featureCount_bin multimapping[True/False]\n" %os.path.realpath(__file__))
		exit()
	
	bam_file = os.path.abspath(argv[1])
	folder = os.path.abspath(argv[2])
	gtf_file = os.path.abspath(argv[3])
	threads = argv[4]
	name = argv[5]
	featureCount_exe = argv[6]
	multimapping= argv[7]

	## Debug
	Debug=True
	
	## variables
	biotype_all(featureCount_exe, folder, gtf_file, bam_file, name, threads, Debug, multimapping)
	## plot results
	RNAbiotypes_stats_file = os.path.join(folder, name + '_RNAbiotype.tsv')
	if files_functions.is_non_zero_file(RNAbiotypes_stats_file):
		pie_plot_results(RNAbiotypes_stats_file, name, folder, Debug)
			
	
######
if __name__== "__main__":
	main()



#===============================================================================
# Example STAR Log Log.final.out
#===============================================================================
#
# 								 Started job on |	Oct 20 10:44:42
# 							 Started mapping on |	Oct 20 10:44:42
# 									Finished on |	Oct 20 11:13:25
# 	   Mapping speed, Million of reads per hour |	63.38
# 
# 						  Number of input reads |	30332431
# 					  Average input read length |	65
# 									UNIQUE READS:
# 				   Uniquely mapped reads number |	11259635
# 						Uniquely mapped reads % |	37.12%
# 						  Average mapped length |	47.60
# 					   Number of splices: Total |	0
# 			Number of splices: Annotated (sjdb) |	0
# 					   Number of splices: GT/AG |	0
# 					   Number of splices: GC/AG |	0
# 					   Number of splices: AT/AC |	0
# 			   Number of splices: Non-canonical |	0
# 					  Mismatch rate per base, % |	0.06%
# 						 Deletion rate per base |	0.00%
# 						Deletion average length |	1.00
# 						Insertion rate per base |	0.00%
# 					   Insertion average length |	1.36
# 							 MULTI-MAPPING READS:
# 		Number of reads mapped to multiple loci |	0
# 			 % of reads mapped to multiple loci |	0.00%
# 		Number of reads mapped to too many loci |	18685152
# 			 % of reads mapped to too many loci |	61.60%
# 								  UNMAPPED READS:
#   Number of reads unmapped: too many mismatches |	268064
# 	   % of reads unmapped: too many mismatches |	0.88%
# 			Number of reads unmapped: too short |	20122
# 				 % of reads unmapped: too short |	0.07%
# 				Number of reads unmapped: other |	99458
# 					 % of reads unmapped: other |	0.33%
# 								  CHIMERIC READS:
# 					   Number of chimeric reads |	0
# 							% of chimeric reads |	0.00%
#===============================================================================
