#!/usr/bin/env python3
############################################################
## Jose F. Sanchez                                        ##
## Copyright (C) 2019-2020 Lauro Sumoy Lab, IGTP, Spain   ##
############################################################
import time
import io
import os
import re
import sys
from io import open
from sys import argv
import pandas as pd
import csv

from HCGB import functions

####################
def generate_DE(dataframe_results, Debug, outfolder):
	"""
	"""
	## get results dictionary for each software employed 
	soft_list = dataframe_results.soft.unique()
	## debugging messages
	if Debug:
		print ("## Debug:")
		print ("soft_list")
		print (soft_list)

	## generate results
	for soft_name in soft_list:
		
		print ("\n+ Summarizing results for software: ", soft_name)
		
		# retrieve files in dictionary
		tmp_df = dataframe_results[ dataframe_results['soft'] == soft_name].set_index('name')
		dict_files = tmp_df['filename'].to_dict()

		if Debug:
			print ("## Debug:")
			print ("tmp_df")
			print (tmp_df)
			print ("dict_files")
			print (dict_files)

		## get data
		(all_data, all_seqs) = generate_matrix(dict_files, soft_name.lower(), Debug)
		
		## discard duplicate UIDs if any
		all_data_filtered, all_data_duplicated = discard_UID_duplicated(all_data)
		
		## dump data in folder provided
		csv_outfile = os.path.join(outfolder, 'miRNA_expression-' + soft_name)
		all_data_filtered.to_csv(csv_outfile + ".csv", quoting=csv.QUOTE_NONNUMERIC)
		all_data_duplicated.to_csv(csv_outfile + '_dup.csv', quoting=csv.QUOTE_NONNUMERIC)
		all_seqs.to_csv(csv_outfile + '_seq.csv', quoting=csv.QUOTE_NONNUMERIC)

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
def generate_matrix(dict_files, soft_name, Debug):
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
		
		## 
		if functions.files_functions.is_non_zero_file(this_file):
			data = pd.read_csv(this_file, sep='\t')
		else:
			print ('\t - Information not available for sample: ', sample)	
			continue
		##
		if (data.size == 0):
			print ('\t - Information not available for sample: ', sample)  
			continue

		## get info, generate unique name and merge for samples
		## header of tsv files: 
		## UID	Read	miRNA	Variant	iso_5p	iso_3p	iso_add3p	iso_snp	sRNAbench

		data['Variant'].fillna('NA', inplace=True)
		data['unique_id'] = data.apply(lambda data: data['miRNA'] + '&' + data['Variant'] + '&' + data['UID'], axis=1)

		## parse according to software
		if (soft_name == 'srnabench'):
			## sRNAbench mirtop creates a column id with sRNAbench instead of sample name
			new_data = data.filter(['unique_id', 'sRNAbench'], axis=1)
			new_data = new_data.set_index('unique_id')
			new_data = new_data.rename(columns={'sRNAbench': sample})

		if (soft_name == 'optimir'):
			## OptimiR mirtop creates a column containing sample name and other tags (trim, joined, fastq...)
			regex=re.compile(sample + '.*')
			search_list = list(filter(regex.match, data.columns.values.tolist()))
			new_data = data.filter(['unique_id', search_list[0]], axis=1)
			new_data = new_data.set_index('unique_id')
			new_data = new_data.rename(columns={search_list[0]: sample})

		if (soft_name == 'miraligner'):
			## miraligner mirtop creates a column containing sample name and other tags (trim, joined, fastq...)
			new_data = data.filter(['unique_id', sample], axis=1)
			new_data = new_data.set_index('unique_id')

		## sequence information
		seq_data = data.filter(['UID', 'Read'], axis=1)	
		seq_data = seq_data.set_index('UID')
		seq_all_data = seq_all_data.append(seq_data, sort=True).drop_duplicates('Read')

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

######
def main():
    ## this code runs when call as a single script

    ## control if options provided or help
    if len(sys.argv) > 1:
        print ("")
    else:
        print ("Provide file containing sample,file")
        exit()        

    dictionary_info = functions.main_functions.file2dictionary(sys.argv[1], ",")
    print (dictionary_info)

    ## get data
    (all_data, all_seqs) = generate_matrix(dictionary_info, "miraligner", False)

    ## discard duplicate UIDs if any
    all_data_filtered, all_data_duplicated = discard_UID_duplicated(all_data)

    ## dump data in folder provided
    outfolder = "./"
    csv_outfile = os.path.join(outfolder, 'miRNA_expression')
    all_data_filtered.to_csv(csv_outfile + ".csv", quoting=csv.QUOTE_NONNUMERIC)
    all_data_duplicated.to_csv(csv_outfile + '_dup.csv', quoting=csv.QUOTE_NONNUMERIC)
    all_seqs.to_csv(csv_outfile + '_seq.csv', quoting=csv.QUOTE_NONNUMERIC)



        
######
if __name__== "__main__":
    main()
    

