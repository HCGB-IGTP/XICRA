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

######
def generate_matrix(dict_files, abs_path_folder):
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
	dict_files = {}
	for l in list_files:
		featurecount_file = abs_path_folder + '/' + l + '/featureCount.out.tsv'
		if os.path.isfile(featurecount_file):
			dict_files[l] = featurecount_file
	
	##	
	all_data = generate_matrix(dict_files, abs_path_folder)

	print ('+ Table contains: ', len(all_data), ' entries\n')
	all_data.to_csv(abs_csv_outfile, quoting=csv.QUOTE_NONNUMERIC)

######
if __name__== "__main__":
	main()
