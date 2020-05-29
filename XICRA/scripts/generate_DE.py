#!/usr/bin/env python3
import time
import io
import os
import re
import sys
from io import open
from sys import argv
import pandas as pd
import csv

from XICRA.scripts import functions
from XICRA.scripts.functions import is_non_zero_file

####################
def generate_DE(dict_files, Debug, outfolder):
	"""
	"""
	## get data
	(all_data, all_seqs) = generate_matrix(dict_files, Debug)
	
	## discard duplicate UIDs if any
	all_data_filtered, all_data_duplicated = discard_UID_duplicated(all_data)
	
	## dump data in folder provided
	csv_outfile = os.path.join(outfolder, 'miRNA_expression')
	all_data_filtered.to_csv(csv_outfile, quoting=csv.QUOTE_NONNUMERIC)
	all_data_duplicated.to_csv(csv_outfile + '_dup', quoting=csv.QUOTE_NONNUMERIC)
	all_seqs.to_csv(csv_outfile + '_seq', quoting=csv.QUOTE_NONNUMERIC)

####################
def discard_UID_duplicated(df_data):
	"""
	"""
	## get data index
	df_data['ID'] = df_data.index
	new_data = df_data.filter(['ID'], axis=1)	

	# split ID (hsa-let-7a-2-3p&NA&qNkjr6Ov2) into miRNA, variant and UID
	tmp = new_data['ID'].str.split('&', expand = True)
	new_data['miRNA']  = tmp[0]
	new_data['variant']  = tmp[1]
	new_data['UID']  = tmp[2]

	## count 
	count_groups = new_data.groupby('UID').count()
	## print to file?
	
	## get duplicated
	bigger1count = count_groups[ count_groups['ID'] > 1 ]
	## print counts to file?
	
	## get list of UIDs duplicate
	bigger1count_list = bigger1count.index.to_list()
	duplicates = new_data[new_data['UID'].isin(bigger1count_list)]
	## print duplicate to file?

	## get duplicated data
	duplicates_indes_list = duplicates.index.to_list()
	duplicates_expression = df_data[df_data.index.isin(duplicates_indes_list)]
	#duplicates_expression['UID'] = duplicates_expression['ID'].str.split('&', expand = True)[2]
	duplicates_expression = duplicates_expression.drop(['ID'], axis=1)
	duplicates_expression.index.name = "ID"
	
	## get clean data
	clean_data_expression = df_data[~df_data.index.isin(duplicates_indes_list)]
	#clean_data_expression['UID'] = clean_data_expression['ID'].str.split('&', expand = True)[2]
	clean_data_expression =	clean_data_expression.drop(['ID'], axis=1)
	clean_data_expression.index.name = "ID"
	## Fix this error
	## A value is trying to be set on a copy of a slice from a DataFrame.
	## Try using .loc[row_indexer,col_indexer] = value instead
	## 
	## See the caveats in the documentation: http://pandas.pydata.org/pandas-docs/stable/indexing.html#indexing-view-versus-copy
	##   clean_data_expression['UID'] = clean_data_expression['ID'].str.split('&', expand = True)[2]
	
	return (clean_data_expression, duplicates_expression)

####################
def generate_matrix(dict_files, Debug):
	"""
	"""
	######################################################
	### For a list of files, generates a count matrix with a unique index id by merging:
	### information provided within each file: name, variant and UID by '&'
	### e.g. AlaAGC&3'-tRF&tRF-16-KSP185D
	### e.g. hsa-let-7a-2-3p&NA&qNkjr6Ov2
	### 
	### It returns a dataframe containing for each index generated
	### count values for each sample (sample_name) in columns
	###
	#########################################################
	all_data = pd.DataFrame()
	seq_all_data = pd.DataFrame()
	for sample, this_file in dict_files.items():
		print ('+ Reading information from sample: ', sample)	
		
		if is_non_zero_file(this_file):
			data = pd.read_csv(this_file, sep='\t')
		else:
			print ('\t - Information not available for sample: ', sample)	
			
		## get info, generate unique name and merge for samples
		## header of tsv files: 
		## UID	Read	miRNA	Variant	iso_5p	iso_3p	iso_add3p	iso_snp	sRNAbench

		data['Variant'].fillna('NA', inplace=True)
		data['unique_id'] = data.apply(lambda data: data['miRNA'] + '&' + data['Variant'] + '&' + data['UID'], axis=1)

		new_data = data.filter(['unique_id', 'sRNAbench'], axis=1)	## change if different from sRNAbench
		new_data = new_data.set_index('unique_id')
	
		seq_data = data.filter(['UID', 'Read'], axis=1)	
		seq_data = seq_data.set_index('UID')
		seq_all_data = seq_all_data.append(seq_data, sort=True).drop_duplicates('Read')

		new_data = new_data.rename(columns={'sRNAbench': sample})
		
		## debugging messages
		if Debug:
			print ("*** DEBUG: data for sample ***")
			print (new_data)
		
		all_data = pd.concat([all_data, new_data], axis=1, sort=True)

	##
	## debugging messages
	if Debug:
		print ("*** DEBUG: data for all samples ***")
		print (all_data)
		print ("*** DEBUG: data for sequences all samples ***")
		print (seq_all_data)
		
		
	return (all_data, seq_all_data)	

######

