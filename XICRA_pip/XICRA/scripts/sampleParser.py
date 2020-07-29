#!/usr/bin/ python3
##########################################################
## Jose F. Sanchez										##
## Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain		##
##########################################################
"""
Prepares samples for further analysis.
"""
## useful imports
import time
import io
import os
import re
import sys
import pandas as pd
from termcolor import colored
import shutil
import concurrent.futures

from XICRA.scripts import functions

###############
def get_fields(file_name_list, pair, Debug, include_all):
	"""
	Get information from files
	
	Creates and returns a dataframe containing information for each sample for:
	"sample", "dirname", "name", "new_name", "name_len", "lane", "read_pair","lane_file","ext","gz", "tag"
	
	:param file_name_list: List of files to parse
	:param pair: True/false for Paired-end files
	:param Debug: True/false for debugging messages
	
	:type file_name_list: list
	:type pair: bool
	:type Debug: bool
	
	:returns: Pandas dataframe.
	
	"""
	## init dataframe
	name_columns = ("sample", "dirname", "name", "new_name", "name_len", 
				"lane", "read_pair","lane_file","ext","gz", "tag", "file")
	name_frame = pd.DataFrame(columns=name_columns)
	
	## loop through list
	for path_files in file_name_list:
		## get file name
		file_name = os.path.basename(path_files)
		dirN = os.path.dirname(path_files)
		trim_search = re.search(r".*trim.*", file_name)
		lane_search = re.search(r".*\_L\d+\_.*", file_name)
		## get name
		if (include_all):
			if (pair):
				name_search = re.search(r"(.*)\_(R1|1|R2|2)\_{0,1}(.*)\.(f.*q)(\..*){0,1}", file_name)
			else:
				name_search = re.search(r"(.*)\.(f.*q)(\..*){0,1}", file_name)

		else:	
			if (pair):
				## pair could be: R1|R2 or 1|2
				if (trim_search):
					name_search = re.search(r"(.*)\_trim\_(R1|1|R2|2)\.(f.*q)(\..*){0,1}", file_name)
					
				else:
					## Lane files: need to merge by file_name: 33i_S5_L004_R1_001.fastq.gz
					## lane should contain L00x			
					if (lane_search):
						name_search = re.search(r"(.*)\_(L\d+)\_(R1|1|R2|2)\_{0,1}(.*)\.(f.*q)(\..*){0,1}", file_name)
					else:
						name_search = re.search(r"(.*)\_(R1|1|R2|2)\.(f.*q)(\..*){0,1}", file_name)
			else:
				if (trim_search):
					name_search = re.search(r"(.*)\_trim(.*)\.(f.*q)(\..*){0,1}", file_name) ## trim.fastq or trim_joined.fastq
				else:
					name_search = re.search(r"(.*)\.(f.*q)(\..*){0,1}", file_name)
		
		### declare
		name= ""
		lane_id= ""
		read_pair= ""
		lane_file= ""
		ext= ""
		gz= ""
	
		if name_search:
			name = name_search.group(1)
			name_len = len(name)
			if (pair):
				if (include_all):
					lane_id = ""
					read_pair = name_search.group(2)
					lane_file = name_search.group(3)
					ext = name_search.group(4)
					gz = name_search.group(5)
								
				## Lane files: need to merge by file_name: 33i_S5_L004_R1_001.fastq.gz
				elif (lane_search):
					lane_id = name_search.group(2)
					read_pair = name_search.group(3)
					lane_file = name_search.group(4)
					ext = name_search.group(5)
					gz = name_search.group(6)
				else:
					## could exist or not
					read_pair = name_search.group(2)
					ext = name_search.group(3)
					gz = name_search.group(4)
			else:
				ext = name_search.group(2)
				gz = name_search.group(3)
	
			## populate dataframe
			name_frame.loc [len(name_frame)] = (path_files, dirN, name, name, name_len, 
											lane_id, read_pair, lane_file, ext, gz, "reads", os.path.basename(path_files))
	
		else:
			## debug message
			if (Debug):
				print (colored("**DEBUG: sampleParser.get_fields **", 'yellow'))
				print (colored("*** ATTENTION: Sample did not match the possible parsing options...", 'yellow'))
				print (file_name)
			
			print (colored("*** ATTENTION: Sample (%s) did not match the possible parsing options..." %path_files, 'yellow'))

	return (name_frame)

###############
def select_samples (list_samples, samples_prefix, pair=True, exclude=False, Debug=False, lane=False, include_all=False):
	"""
	Select samples
	
	Given a sample prefix (any or a given list), this function retrieves
	sample files from a list given. If exclude option provided, excludes 
	the files retrieved from the total.
	
	:param list_samples: List of absolute path for fastq files
	:param samples_prefix: List of prefix to search 
	:param pair: True/false for paired-end files
	:param exclude: True/false for exclude found files from the total
	:param Debug: True/false for debugging messages
	:param lane: Include lane tag within name id
	
	:type list_samples: list
	:type samples_prefix: list
	:type pair: bool
	:type exclude: bool
	:type Debug: bool
	:type lane: bool
	
	:returns: Dataframe
	"""
	#Get all files in the folder "path_to_samples"
	sample_list = pd.DataFrame(columns=('sample', 'file'))
	
	for names in samples_prefix:
		for path_fastq in list_samples:	
			fastq = os.path.basename(path_fastq)
			samplename_search = re.search(r"(%s)\_{0,1}(R1|1|R2|2){0,1}(.*){0,1}\.f.*q.*" % names, fastq)
			enter = ""
			if samplename_search:
				if (exclude): ## exclude==True
					enter = False
				else: ## exclude==True
					enter = True
			else:
				if (exclude): ## exclude==True
					enter = True
				else: ## exclude==True
					enter = False
					
			if enter:
				if fastq.endswith('.gz') or fastq.endswith('fastq') or fastq.endswith('fq'):
					sample_list.loc[len(sample_list)] = (names, path_fastq) 
				else:
					## debug message
					if (Debug):
						print (colored("**DEBUG: sampleParser.select_samples **", 'yellow'))
						print (colored("** ERROR: %s is a file that is neither in fastq.gz or .fastq format, so it is not included" %path_fastq, 'yellow'))
							
	## discard duplicates if any
	non_duplicate_names = sample_list['sample'].to_list() #
	non_duplicate_names = list(set(non_duplicate_names))
	
	## it might be a bug in exclude list.
	## if sample X1 is provided to be excluded, we might be also excluding
	## sample X12, sample X13, etc.
	## TODO: check this

	## debugging messages
	if Debug:
		print (colored("** DEBUG: select_samples",'yellow'))
		print ("non_duplicate_names:")
		print (non_duplicate_names)
	
	## check they match with given input
	if (exclude): ## exclude==True
		if bool(set(samples_prefix).intersection(non_duplicate_names)):
			print(colored("** ERROR: Some non desired samples are included", 'red'))
	else: ## exclude==True
		non_duplicate_names = set(samples_prefix).intersection(non_duplicate_names)

	## get fields
	
	tmp = sample_list[ sample_list['sample'].isin(non_duplicate_names) ]
	non_duplicate_samples = tmp['file'].to_list()
	
	## debugging messages
	if Debug:
		print (colored("** DEBUG: select_samples",'yellow'))
		print ("non_duplicate_names:")
		print (non_duplicate_names)
		print ("samples_prefix")
		print (samples_prefix)
		print ("non_duplicate_samples")
		print (non_duplicate_samples)
		print ("tmp dataframe")
		#functions.print_all_pandaDF(tmp)
		print(tmp)
				
	## get info
	name_frame_samples = get_fields(non_duplicate_samples, pair, Debug, include_all)	
	number_files = name_frame_samples.index.size
	total_samples = set(name_frame_samples['name'].to_list())
	
	##
	if (lane):
		## include lane tag within name
		name_frame_samples['name'] = name_frame_samples['name'] + '_' + name_frame_samples['lane']
		name_frame_samples['new_name'] = name_frame_samples['name']
			
	## debugging messages
	if Debug:
		print (colored("** DEBUG: select_samples",'yellow'))
		print ("name_frame_samples:")
		print (name_frame_samples)
		print ("number_files:")
		print (number_files)
		print ("total_samples:")
		print (total_samples)
	
	### get some stats
	if (number_files == 0):
		print (colored("\n**ERROR: No samples were retrieved. Check the input provided\n",'red'))
		exit()
	print (colored("\t" + str(number_files) + " files selected...", 'yellow'))
	print (colored("\t" + str(len(total_samples)) + " samples selected...", 'yellow'))
	if (pair):
		print (colored("\tPaired-end mode selected...", 'yellow'))
	else:
		print (colored("\tSingle end mode selected...", 'yellow'))
	
	## return info
	return (name_frame_samples)

###############
	

###############    
def gunzip_merge(outfile, list_files):
	"""
	Merge gunzip files into final file
	
	:param outfile: String for output file
	:param list_files: List of files to merge
	
	:type outfile: string
	:type list_files: list
		
	"""
	list_files = list(list_files)
	list_files.sort()
	print ("\tMerging files into: ", outfile)
	#print ("\tFiles: ", list_files)

	with open(outfile, 'wb') as wfp:
		for fn in list_files:
			with open(fn, 'rb') as rfp:
				shutil.copyfileobj(rfp, wfp)

	return()
	
###############    
def one_file_per_sample(dataFrame, outdir_dict, threads, outdir, Debug=False):
	"""
	Merge fastq files from different lanes positions for each sample
	
	"""
	## merge sequencing files for sample, no matter of sector or lane generated.	
	list_samples = set(dataFrame['new_name'].tolist())
	print (colored("\t" + str(len(list_samples)) + " samples to be merged from the input provided...", 'yellow'))
	print ("+ Merging sequencing files for samples")

	##
	sample_frame = dataFrame.groupby(["new_name", "read_pair"])
	
	### get extension for files
	ext_list = dataFrame.ext.unique()
	gz_list = dataFrame.gz.unique()
	## might generate a bug if several extension or some zip/unzip files provided
	if gz_list:
		ext = ext_list[0] + gz_list[0]
	else:
		ext = ext_list[0]

	# We can use a with statement to ensure threads are cleaned up promptly
	with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor: ## need to do 1 by one as there is a problem with the working directory
		commandsSent = { executor.submit(gunzip_merge, 
										outdir_dict[name[0]] + '/' + name[0] + '_' + name[1] + '.' + ext, 
										sorted(set(cluster["sample"].tolist()))): name for name, cluster in sample_frame }
		for cmd2 in concurrent.futures.as_completed(commandsSent):
			details = commandsSent[cmd2]
			try:
				data = cmd2.result()
			except Exception as exc:
				print ('***ERROR:')
				print (cmd2)
				print('%r generated an exception: %s' % (details, exc))				
							
	## return output name merged generated in dataframe
	name_columns = ("new_name", "dirname", "read_pair", "new_file", "ext", "gz")
	name_frame = pd.DataFrame(columns=name_columns)

	## print to a file
	timestamp = functions.create_human_timestamp()
	merge_details = outdir + '/' + timestamp + '_prep_mergeDetails.txt'
	merge_details_hd = open(merge_details, 'w')

	for name, cluster in sample_frame: ## loop over samples
		outfile = outdir_dict[name[0]] + '/' + name[0] + '_' + name[1] + ext
		
		merge_details_hd.write("####################\n")		
		merge_details_hd.write("Sample: " + name[0] + '\n')
		merge_details_hd.write("New name: " + name[0] + '\n')
		
		merge_details_hd.write("Read: " + name[1] + '\n')
		merge_details_hd.write("Files:\n")
		merge_details_hd.write(",".join(sorted(cluster["sample"].tolist())))
		merge_details_hd.write('\n')
		merge_details_hd.write("####################\n")		
		
		name_frame.loc[len(name_frame)] = (name[0], outdir_dict[name[0]], name[1], outfile, ext_list[0], gz_list[0])

	merge_details_hd.close()
	return(name_frame)

################################
def get_files(options, input_dir, mode, extension):
	"""
	Parser function to get files
	
	Given an input dir and a mode in retrieves
	matching files with the extension desired.
	
	:param options: XICRA options as parser.parse_args options.
	:param input_dir: Absolute path to input dir containing samples.
	:param mode: Options are: fastq, trim, annot, assembly.
	:param extension: List of possible extension to retrieve.
	
	:type options: parser.parse_args
	:type input_dir: string 
	:type mode: string 
	:type extension: list
	
	:returns: Pandas dataframe with sample and file information.
	"""
	## get list of input files
	files = []
	print ()
	functions.print_sepLine("-",50, False)
	print ('+ Getting files from input folder... ')
	print ('+ Mode: ', mode,'. Extension:', extension)
	if (options.project):
		
		## input folder is not a dir, is it a batch input file?
		if (options.batch):
			if os.path.isfile(input_dir):
				if (options.debug):
					print (colored("\n**DEBUG: sampleParser.get_files input folder is a batch file, get full path **", 'yellow'))
				dir_list = [line.rstrip('\n') for line in open(input_dir)]
				for d in dir_list:
					if os.path.exists(d):
						print ('+ Folder (%s) exists' %d)
						files = files + functions.get_fullpath_list(d)
					else:
						## input folder does not exist...
						if (options.debug):
							print (colored("\n**DEBUG: sampleParser.get_files batch option; input folder does not exists **", 'yellow'))
							print (d)
							print ("\n")
		else:		
			### a folder containing a project is provided
			if os.path.exists(input_dir):
				print ('+ Input folder exists')
				## get files in folder
				for ext in extension:
					if mode == 'trim':
						files_tmp = functions.get_fullpath_list(input_dir)
						files = [s for s in files_tmp if ext in s]
					else:
						files_tmp = functions.retrieve_matching_files(input_dir, ext)				
						files = files + files_tmp
	
				files = set(files)

			else:
				## input folder does not exist...
				if (options.debug):
					print (colored("\n**DEBUG: sampleParser.get_files input folder does not exists **", 'yellow'))
					print (input_dir)
					print ("\n")
	
				print (colored('***ERROR: input folder does not exist or it is not readable', 'red'))
				exit()						
	else:
		### provide a single folder or a file with multiple paths (option batch)
		if (options.batch):
			if os.path.isfile(input_dir):
				dir_list = [line.rstrip('\n') for line in open(input_dir)]
				for d in dir_list:
					if os.path.exists(d):
						print ('+ Folder (%s) exists' %d)
						files_tmp = functions.get_fullpath_list(d)
						files = files + files_tmp
			
					else:
						## input folder does not exist...
						if (options.debug):
							print (colored("\n**DEBUG: sampleParser.get_files batch option; input folder does not exists **", 'yellow'))
							print (d)
							print ("\n")
						
		else:
			if os.path.exists(input_dir):
				print ('+ Input folder exists')
				## get files in folder
				files = functions.get_fullpath_list(input_dir)
			else:
				## input folder does not exist...
				if (options.debug):
					print (colored("\n**DEBUG: sampleParser.get_files input folder does not exists **", 'yellow'))
					print (input_dir)
					print ("\n")
	
				print (colored('***ERROR: input folder does not exist or it is not readable', 'red'))
				exit()

	if options.debug:
		print (colored("** DEBUG: sampleParser.get_files files", 'yellow'))
		print (files)

	## get list of samples
	samples_names = []
	exclude=False

	if (options.in_sample):
		if os.path.isfile(os.path.abspath(options.in_sample)):
			in_file = os.path.abspath(options.in_sample)
			samples_names = [line.rstrip('\n') for line in open(in_file)]
			print ('+ Retrieve selected samples to obtain from the list files available.')
		else:
			sample_names = options.in_sample
		
		exclude=False

		## in sample list...
		if (options.debug):
			print (colored("\n**DEBUG: sampleParser.get_files include sample list **", 'yellow'))
			print (samples_names, '\n')


	elif (options.ex_sample):
		ex_file = os.path.abspath(options.ex_sample)
		samples_names = [line.rstrip('\n') for line in open(ex_file)]
		print ('+ Retrieve selected samples to exclude from the list files available.')		
		exclude=True

		## in sample list...
		if (options.debug):
			print (colored("\n**DEBUG: sampleParser.get_files exclude sample list **", 'yellow'))
			print (samples_names, '\n')

	else:
		samples_names = ['.*']

	## discard empty sample_names
	samples_names = list(filter(None, samples_names)) ## empty space

	## discard some files obtain
	files = [s for s in files if '.bam' not in s]
	files = [s for s in files if '.sam' not in s]
	files = [s for s in files if '.log' not in s]
	files = [s for s in files if '.annot' not in s]
	files = [s for s in files if '.abundances.txt' not in s]
	files = [s for s in files if '.gff3' not in s]
	files = [s for s in files if 'trimmed.fq' not in s]
	files = [s for s in files if 'trim.clpsd.fq' not in s]
	files = [s for s in files if 'failed.fq.gz' not in s]
	files = [s for s in files if 'unjoin' not in s]
	
	if (mode == 'fastq'):
		files = [s for s in files if 'trim' not in s]
	
	##
	files = list(filter(None, files)) ## empty space
		
	## files list...
	if (options.debug):
		print (colored("\n**DEBUG: sampleParser.get_files files list to check **", 'yellow'))
		print ('DO NOT PRINT THIS LIST: It could be very large...')
		##print (files, '\n')

	## get information
	if mode in ['fastq', 'trim', 'join']:
		pd_samples_retrieved = select_samples(files, samples_names, options.pair, exclude, options.debug, options.include_lane, options.include_all)
	else:
		pd_samples_retrieved = select_other_samples(options.project, files, samples_names, mode, extension, exclude, options.debug)		
		
	return(pd_samples_retrieved)

