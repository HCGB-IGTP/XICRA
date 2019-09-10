#usr/bin/env python
import time
import io
import os
import re
import sys
from io import open
from sys import argv
import pandas as pd
import csv

######
def generate_matrix(list_files, abs_path_folder):
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
	for f in list_files:
		this_file = abs_path_folder + '/' + f

		print ('+ Reading information from file: ', this_file)	
		data = pd.read_csv(this_file, sep='\t')
	
		## skip if file is empty
		if data.empty:
			continue

		## get info, generate unique name and merge for samples
		data['variant'].fillna('NA', inplace=True)
		data['unique_id'] = data.apply(lambda data: data['name'] + '&' + data['variant'] + '&' + data['UID'], axis=1)

		new_data = data.filter(['unique_id', 'expression'], axis=1)	
		new_data = new_data.set_index('unique_id')
	
		sample_name = data['sample_name'].to_list()
		new_data = new_data.rename(columns={'expression': sample_name[0]})
		all_data = pd.concat([all_data, new_data], axis=1, sort=True)

	##
	return (all_data)	

######
def main():
	## ARGV
	if len (sys.argv) < 3:
		print ("\nUsage:")
		print ("python3 %s folder out_csv\n" %os.path.realpath(__file__))
		exit()

	abs_path_folder = os.path.abspath(argv[1])
	abs_csv_outfile = os.path.abspath(argv[2])
	list_files = os.listdir(abs_path_folder)
	all_data = generate_matrix(list_files, abs_path_folder)

	print ('+ Database contains: ', len(all_data), ' entries\n')
	all_data.to_csv(abs_csv_outfile, quoting=csv.QUOTE_NONNUMERIC)

######
if __name__== "__main__":
	main()
