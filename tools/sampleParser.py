#usr/bin/env python

## useful imports
import time
import io
import os
import re
import subprocess
import sys
from sys import argv

## import functions
toolDir = os.path.dirname(argv[0])
sys.path.append(toolDir)
import functions

###############
def select_samples (samples_prefix, path_to_samples):
    
    #Get all files in the folder "path_to_samples"    
	files = os.listdir(path_to_samples)
	sample_list = []
	for fastq in files:	
		samplename_search = re.search(r"(%s)\_.*" % samples_prefix, fastq)
		if samplename_search:
			if 'merged' not in fastq:
				if fastq.endswith('.gz'):
					sample_list.append(fastq)
				elif fastq.endswith('fastq'):
					sample_list.append(fastq)
				else:
					print ("** ERROR: ", fastq, 'is a file that is neither in fastq.gz or .fastq format, so it is not included')

	non_duplicate_samples = list(set(sample_list))	
	discard_samples = []
	for files in non_duplicate_samples:
		if files.endswith('.gz'):
			gz_search = re.search(r"(.*)\.gz", files)
			if gz_search:
				file_name = gz_search.group(1)
				if file_name not in non_duplicate_samples:
					functions.extract(path_to_samples + '/' + files)
					non_duplicate_samples.append(file_name)

			discard_samples.append(files)			
	
	for samples in discard_samples:
		non_duplicate_samples.remove(samples)
						
	non_duplicate_samples2 = list(set(non_duplicate_samples))						
	number_samples = len(non_duplicate_samples2)	

	print ("\t\t- ", number_samples," samples selected from ", path_to_samples)
	return sorted(non_duplicate_samples)
	
###############

###############    
def one_file_per_sample(final_sample_list, path_to_samples, directory, read, output_file, prefix_list, num_threads):
	## merge sequencing files for sample, no matter of sector or lane generated.
	grouped_subsamples = []
	bigfile_list = []
	commands2sent = []
	
	output_file = open(output_file, 'a')
	output_file.write("\nMerge samples:\n")
	
	for samplex in final_sample_list:
		if samplex not in grouped_subsamples:
			samplename = []
			subsamples = []
			for prefix in prefix_list:
				samplename_search = re.search(r"(.*)\_(%s)" % read, samplex)
				if samplename_search:			
					original_name = samplename_search.group(1)
					name_search = re.search(r".*%s.*" % prefix, original_name)
					if name_search:
						commonname = original_name + "_" + read + ".fastq"
						bigfilepath = directory + "/" + commonname
						bigfile_list.append(commonname)	
					
						for sampley in final_sample_list:
							if original_name in sampley:
								subsamples.append(path_to_samples + "/" + sampley)
								grouped_subsamples.append(sampley)
						if not os.path.isfile(bigfilepath) or os.stat(bigfilepath).st_size == 0:
							partsofsample = ' '.join(sorted(subsamples))
							cmd = 'cat %s >> %s' %(partsofsample, bigfilepath)
							## DUMP in file					
							output_file.write(cmd)   
							output_file.write('\n')						
							## get command				
							commands2sent.append(cmd)
						else:
							print ('\t + Sample %s is already merged' % commonname)
	## close file
	output_file.close()		
	
	#sent commands on threads
	functions.sender(commands2sent, num_threads)	
	print ('There are' , len(bigfile_list) , 'samples after merging for read' , read, '\n')
	return bigfile_list
###############

