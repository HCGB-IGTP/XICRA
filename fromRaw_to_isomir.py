#usr/bin/env python

## useful imports
import sys
from sys import argv
import time
from datetime import datetime
from io import open
import subprocess
import os
import re

from configparser import ConfigParser
#import concurrent.futures

#####################
#### functions ######
#####################

def help_options():
	print ("\n#######################################################################")
	print ("  NAME: fromRaw_to_isomiR")
	print ("  VERSION: 0.3")
	print ("  AUTHORS: Antonio Luna de Haro (v0.1) & Jose F Sanchez-Herrero (v1).")
	print ("           Copyright (C) 2018-2019 Lauro Sumoy Lab, IGTP, Spain")
	print ("#########################################################################")
	print ("\nDESCRIPTION:")
	print ("- This script is a pipeline generated for the trimming and joining of paired-end or single end reads from small RNA-seq data.")
	print ("\n- Available analysis are:")
	print ("\t+ miRNA-isomiR analysis using sRNAtoolbox ")
	print ("\t+ tRFs using MINTmap and MINTbase.")
	print ("")
	print ("USAGE:\npython", os.path.abspath(argv[0]),"config_file.txt ")
	print ("\nPARAMETERS:")
	print ("A configuration file is necessary that includes general and detailed information for the project.")
	print ("")
	print ("*******************************************************")
	print (" Configuration file details:")
	print ("*******************************************************")
	print ("[GENERAL]")
	print ("fastq_R1 = /path/to/file/fastqR1")
	print ("fastq_R2 = /path/to/file/fastqR2")
	print ("project = project_name")
	print ("")
	print ("[VARIABLES]")
	print ("prefix = prefix_name_sample_selection")
	print ("thread = num_threads")
	print ("option = isomiR|tRFs|all")
	print ("merge_samples = YES|NO")
	print ("")
	print ("[PARAMETERS]")
	print ("adapter_3 = sequence1")
	print ("adapter_5 = sequence2")
	print ("fastqjoin_percent_difference = 8")
	print ("")
	print ("[EXECUTABLES]")
	print ("cutadapt = /path/to/cutadapt/bin/cutadapt")
	print ("fastqjoin = /path/to/fastqjoin_path/fastq-join")
	print ("sRNAbenchtoolbox = /path/to/sRNAtoolboxDB_folder")
	print ("mirtop_exec /path/to/mirtop_bin/mirtop")
	print ("MINTmap_folder = /path/to/MINTmap/folder")
	print ("")
	print ("[miRNA_ANALYSIS]")
	print ("miRNA_gtf = /path/to/human_genome/gff_file/for/miRNA/file.gff3")
	print ("")
	print ("*******************************************************")
	print ("")
	print ("")
	print ("CITATION:")
	print ("[to add citation]")
	print ("")
	print ("")
	print ("DOCUMENTATION:")
	print ("See [ <http://website/> ] for full documentation")
	print ("")
	print ("")
	print ("LICENSE:")
	print ("This program is free software: you can redistribute it and/or modify")
	print ("it under the terms of the GNU General Public License as published by")
	print ("the Free Software Foundation, either version 3 of the License, or")
	print ("(at your option) any later version.")
	print ("")
	print ("This program is distributed in the hope that it will be useful,")
	print ("but WITHOUT ANY WARRANTY; without even the implied warranty of")
	print ("MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the")
	print ("GNU General Public License for more details.")
	print ("")
	print ("You should have received a copy of the GNU General Public License")
	print ("along with this program.  If not, see <http://www.gnu.org/licenses/>.")
	print ("")
	print ("")
	print ("#################################################\n\n")
###############

###############   
def getime (start_time):
    total_sec = time.time() - start_time
    m, s = divmod(int(total_sec), 60)
    h, m = divmod(m, 60)
    return h, m, s
###############   
    
###############   
def create_subfolder (name, path):
    ## create subfolder  ##	
	subfolder_path = path + "/" + name
	    
    # define the access rights
	try:
		os.mkdir(subfolder_path, access_rights)
	except OSError:  
	   	print ("\tDirectory %s already exists" % subfolder_path)
	else:  
		print ("\tSuccessfully created the directory %s " % subfolder_path)
	
	print ("")
	return subfolder_path
###############   
    
###############   
def timestamp (start_time_partial):
	h,m,s = getime(start_time_partial)
	print ('--------------------------------')
	print ('(Time spent: %i h %i min %i s)' %(int(h), int(m), int(s)))
	print ('--------------------------------')
	return time.time()
############### 

###############
def select_samples (samples_prefix, path_to_samples):
    
    #Get all files in the folder "path_to_samples"    
	files = os.listdir(path_to_samples)
	sample_list = []
	if samples_prefix == 'all':
		for fastq in files:
			if 'merged' in fastq: 
				continue
			
			## ToDOs: check file is unzipped
			if fastq.endswith('.gz'):
				sample_list.append(fastq[:-3]) 
			elif fastq.endswith('fastq'):
				sample_list.append(fastq)
			else:
				print ("** ERROR: ", fastq, 'is a file that is neither in fastq.gz or .fastq format, so it is not included')
	else:
		for fastq in files:
			if fastq.startswith(samples_prefix) and 'merged' not in fastq:
				if fastq.endswith('.gz'):
					sample_list.append(fastq[:-3])
				elif fastq.endswith('fastq'):
					sample_list.append(fastq)
				else:
					print ("** ERROR: ", fastq, 'is a file that is neither in fastq.gz or .fastq format, so it is not included')

									
	non_duplicate_samples = list(set(sample_list))
	number_samples = len(non_duplicate_samples)
	print ("\t\t- ", number_samples," samples selected from ", path_to_samples)
	return sorted(non_duplicate_samples)
###############
 
###############    
def one_file_per_sample(final_sample_list, path_to_samples, directory, read):
	## merge sequencing files for sample, no matter of sector or lane generated.
	grouped_subsamples = []
	bigfile_list = []

	for samplex in final_sample_list:
		if samplex not in grouped_subsamples:
			#print "\nTest: ", samplex
			subsamples = []
			samplename_search = re.search('([a-zA-Z]{2,3})\_(\d{1,2})\_([a-zA-Z]{6})(.*)', samplex)
			if samplename_search:
				path = samplename_search.group(1)
				sample = samplename_search.group(2)
				comonpart = path + "_" + sample + "_"
				commonname = path + "_" + sample + "_" + read + ".fastq"
				bigfilepath = directory + "/" + commonname
				bigfile_list.append(commonname)			
				
				for sampley in final_sample_list:
					if comonpart in sampley:
						subsamples.append(path_to_samples + "/" + sampley)
						grouped_subsamples.append(sampley)
				if not os.path.isfile(bigfilepath) or os.stat(bigfilepath).st_size == 0:
					#partsofsample = ' '.join(sorted(subsamples))
					partsofsample=subsamples[0]
					cmd = 'cat %s >> %s' %(partsofsample, bigfilepath)
					
					## ToDOs: DUMP in file merge_info.txt
					try:
						print ('\t+ %s created' %commonname)
						subprocess.check_output(cmd, shell = True)
					except subprocess.CalledProcessError as err:
						print (err.output)
						sys.exit()
				else:
					print ('\t + Sample %s is already merged' % commonname)
	print ('There are' , len(bigfile_list) , 'samples after merging for read' , read, '\n')
	return bigfile_list
	
###############

###############
def get_symbolic_link (final_sample_list, path_to_samples, directory):
	for samplex in final_sample_list:
		sample_path = path_to_samples + '/' + samplex
		cmd = 'ln -s %s %s' %(sample_path, directory)
		try:
			print ('Symbolic link sample ', samplex)
			subprocess.check_output(cmd, shell = True)
		except subprocess.CalledProcessError as err:
			print (err.output)
			#sys.exit()

	files2return = os.listdir(directory)
	return files2return
###############

###############   
def cutadapt (list_R1, list_R2, path, out_path):
	trimmed_R1 = []
	trimmed_R2 = [] 

	for file_R1 in list_R1:
		file_R1_path = path + "/" + file_R1
		cmd = []
		cutadapt_exe = config['EXECUTABLES']['cutadapt']    	
		adapter_3 = config['PARAMETERS']['adapter_3']
		
		sampleR1_search = re.search('([a-zA-Z]{2,3})\_(\d{1,2})\_.*', file_R1)	
		if sampleR1_search:
			name = sampleR1_search.group(1) + '_' + sampleR1_search.group(2)
			path_name = out_path + "/" + name
			common = path_name + '_trimmed_'
			o_param = common + "R1.fastq"
			p_param = common + "R2.fastq"
			trimmed_R1.append(o_param)

			if paired_end == False:
				## single-end
				cmd = '%s -a %s -o %s %s' %(cutadapt_exe, adapter_3, o_param, file_R1_path)
			else: ## paired-end
				if list_R2: ## paired-end
					sampleR2 = path + "/" + sampleR1_search.group(1) + '_' + sampleR1_search.group(2) + '_R2.fastq'
					try:
						os.path.isfile(sampleR2)					
					except:
						print ("**ERROR: pair for sample ",o_param," doest not exist")
						print ("Sample will be treated as single end")
						## set cmd for single eng as no R2 file
						cmd = '%s -a %s -o %s %s' %(cutadapt_exe, adapter_3, o_param, file_R1_path)
					else:
						#paired end:
						adapter_5 = config['PARAMETERS']['adapter_5']
						file_R2_path = sampleR2
						cmd = '%s -a %s -A %s -o %s -p %s %s %s' %(cutadapt_exe, adapter_3, adapter_5, o_param, p_param, file_R1_path, file_R2_path)
						trimmed_R2.append(p_param)
				else:
					## single-end
					cmd = '%s -a %s -o %s %s' %(cutadapt_exe, adapter_3, o_param, file_R1_path)

			if (os.path.isfile(o_param)):
				if (os.path.isfile(p_param)):
					print ('\tSample %s is already trimmed in paired-end mode' % name)
				else:
					print ('\tSample %s is already trimmed in single-end mode' % name)		

			## not trimmed
			else: 
				## ToDOs: DUMP in file cutadapt_info.txt
				#print 'The following cmd is being executed at the shell: \n', cmd			
				try:
					#print 'Trimming for sample ', name
					subprocess.check_output(cmd, shell = True)
					print ('\tAdapters trimmed for the sample ', name)
				except subprocess.CalledProcessError as err:
					print (err.output)
					sys.exit()

	return trimmed_R1, trimmed_R2
###############

###############     
def fastqjoin (trimmed_R1, trimmed_R2, out_path):
	joined_files = []
	cmd = []
	fastqjoin_exe = config['EXECUTABLES']['fastqjoin']    	
	error_param = config['PARAMETERS']['fastqjoin_percent_difference']
	for file_R1 in trimmed_R1:
		sampleR1_search = re.search('([a-zA-Z]{2,3})\_(\d{1,2})\_trimmed_.*', file_R1)	
		if sampleR1_search:
			# names
			name = sampleR1_search.group(1) + '_' + sampleR1_search.group(2)
			path_name = out_path + '/' + name + '_trimmed_join.fastq'
			unjoined_1 = out_path + '/' + name + '_trimmed_unjoin_R1.fastq'
			unjoined_2 = out_path + '/' + name + '_trimmed_unjoin_R2.fastq'
			joined_files.append(path_name)
			
			## get read2
			string2search = '.*' + name + '.*'
			regex=re.compile(string2search)			
			sampleR2 = [m.group(0) for l in trimmed_R2 for m in [regex.search(l)] if m]

			##
			if sampleR2[0] in trimmed_R2:
				cmd = fastqjoin_exe + ' -p %s %s %s -o %s -o %s -o %s' %(error_param, file_R1, sampleR2[0], unjoined_1, unjoined_2, path_name)
			else: ## single-end
				cmd = 'cp %s %s' %(file_R1, path_name)
			if (os.path.isfile(path_name)):
				print ('\tSample %s is already joined' %name)
			else: ## not merged
				## ToDOs: print in file fastqjoin_info.txt
				#print 'The following cmd is being executed at the shell: \n', cmd
				try:
					print ('\tMerge sample ', name)
					subprocess.check_output(cmd, shell = True)
				except subprocess.CalledProcessError as err:
					print (err.output)
					sys.exit()
			
			## ToDOs: count and provide statistics for joined reads
					
	return joined_files
###############       
    
###############       
def sRNAbench (joined_reads, outpath):
	results = []
	sRNAbench_exe = config['EXECUTABLES']['sRNAbenchtoolbox'] + 'exec/sRNAbench.jar'
	sRNAbench_db = config['EXECUTABLES']['sRNAbenchtoolbox'] 	
	for jread in joined_reads:
		sample_search = re.search('.*\/([a-zA-Z]{2,3})\_(\d{1,2})\_trimmed.*', jread)	
		if sample_search:
			outdir = sample_search.group(1) + "_" + sample_search.group(2)
			finalpath = outpath + '/' + outdir + '/'
			results.append(finalpath)
			if (os.path.isdir(finalpath)):
				print ('\tisomiRs analysis for sample %s already exists' %outdir)
			else:
				cmd = 'java -jar %s dbPath=%s input=%s output=%s microRNA=hsa isoMiR=true plotLibs=true graphics=true plotMiR=true bedGraphMode=true writeGenomeDist=true chromosomeLevel=true chrMappingByLength=true ' %(sRNAbench_exe, sRNAbench_db, jread, finalpath)
				# ToDOs: print into file sRNAbench
				#print 'The following cmd is being executed at the shell: \n', cmd
				try:
					subprocess.check_output(cmd, shell = True)
					print ('\tChecked sample %s for isomiRs' %outdir)

				except subprocess.CalledProcessError as err:
					print (err.output)
					sys.exit()
	return results
###############   
    
###############   
def miRTop (results, outpath):
	gtfs = []
	sRNAbench_hairpin = config['EXECUTABLES']['sRNAbenchtoolbox'] + 'libs/hairpin.fa'
	mirtop_exec = config['EXECUTABLES']['mirtop_exec']
	miRNA_gtf = config['miRNA_ANALYSIS']['miRNA_gtf']	
	species = 'hsa' #homo sapiens
	for folder in results:
		name = folder.split('/')[-2]
		outdir = outpath + '/' + name
		outdir_stats = outdir + "/Stats"
		outdir_gtf = outdir + "/mirtop.gtf"
		if (os.path.isdir(outdir)):
			print ('\tSample %s has already an isomiRs gtf file' %name)
		else:
			cmd = mirtop_exec + ' gff --sps %s --hairpin %s --gtf %s --format srnabench -o %s %s' %(species, sRNAbench_hairpin, miRNA_gtf, outdir, folder)

			try:
				print ("\n##########################################")
				print ('Creating isomiRs gtf file for sample %s' %name)
				print ("-------------------------------------------")
				# ToDOs: print to file mirtop_info.txt
				print ('The following cmd is being executed at the shell: \n', cmd)
				print ("")
				subprocess.check_output(cmd, shell = True)

			except subprocess.CalledProcessError as err:
				print (err.output)
				sys.exit()

		## get stats for each
		if (os.path.isdir(outdir_stats)):
			print ('\tSample %s has already an isomiRs stats folder' %name)
		else:
			cmd_stats = mirtop_exec + ' stats -o %s %s' %(outdir_stats, outdir_gtf)
			try:
				print ("\n##########################################")
				print ('Creating isomiRs stats for sample %s' %name)
				print ("-------------------------------------------")
				# ToDOs: print to file mirtop_info.txt
				print ('The following cmd is being executed at the shell: \n', cmd_stats)
				print ("##########################################")
				subprocess.check_output(cmd_stats, shell = True)
			except subprocess.CalledProcessError as err:
				print (err.output)
				#sys.exit()
    
###############

###############
def isomiR_analysis (path, count, reads, time_partial):
	##############################
    ####### Step: sRNAbench ######
	##############################
	print ("\n+ Run sRNAbenchtoolbox:")
	name_sRNAbench_folder = str(count) + '.isomiR_sRNAbenchtoolbox'
	sRNAbench_folder = create_subfolder(name_sRNAbench_folder, path)
	results = sRNAbench(reads, sRNAbench_folder)
	
	## timestamp
	time_partial = timestamp(time_partial)
	
	############################
    ####### Step: miRTop #######
	############################
	name_miRTop_folder = str(count) + '.isomiR_miRTop'
	miRTop_folder = create_subfolder(name_miRTop_folder, path)
	gtfs=miRTop(results, miRTop_folder)
	
	## timestamp
	time_partial = timestamp(time_partial)

	##############################################
    ####### Step: create expression matrix #######
	##############################################
	name_isomiR_matrix_folder = str(count) + '.isomiR_matrix'
	isomiR_matrix_folder = create_subfolder(name_isomiR_matrix_folder, path)

	for gtffile in gtfs:

		sample = gtffile.rpartition('/')[-1][:-4]
		print ("\tParsing sample ", sample)

		## parse gtf file
		sample_dict = parse_gtf(gtffile)

		## get filename
		filename = isomiR_matrix_folder + '/' + sample
		filenames.append(filename)

		## create matrix
		make_table(sample_dict, filename)
		print ('\tExpression matrix created for ', sample)
	
	## timestamp
	time_partial = timestamp(time_partial)

	return time_partial
	
###############

###############   
def parse_gtf(gtffile):
    gtfile = open(gtffile)
    text = gtfile.read()
    lines = text.splitlines()
    sample_dict = {}
    
    for line in lines:
        if not line.startswith('#'):
            name = line.split('\t')[0]
            name = line.split('\t')[-1].split(';')[2].split()[-1]
            variant = line.split('\t')[-1].split(';')[4].split()[-1]
            variant = '-' + variant
            if variant == '-NA':
                variant = ''
            expression = int(line.split('\t')[-1].split(';')[6].split()[-1])
            namevariant = name + variant
            if namevariant in sample_dict.keys():
                old = sample_dict[namevariant]
                new = old + expression
                sample_dict[namevariant] = new
            else:
                sample_dict[namevariant] = expression
            #sample_dict[name+variant]= {}
    return sample_dict
###############

###############     
def make_table(dictionary, filename):
    fil = open(filename, 'w')    
    isomirs = dictionary.keys()
    
    fil.write('Isomir\t')
    fil.write(filename.partition('_expression')[0])
    fil.write('\n')
    for sample in isomirs:
        fil.write(sample)
        fil.write('\t')
        fil.write(str(dictionary[sample]))
        fil.write('\n')
    fil.close()      
###############

###############
def MINTmap(reads, folder):
	MINTmap = config['EXECUTABLES']['MINTmap_folder'] + 'MINTmap.pl'	
	MINTmap_table = config['EXECUTABLES']['MINTmap_folder'] + 'LookupTable.tRFs.MINTmap_v1.txt'
	MINTmap_tRNAseq = config['EXECUTABLES']['MINTmap_folder'] + 'tRNAspace.Spliced.Sequences.MINTmap_v1.fa'
	MINTmap_tRF = config['EXECUTABLES']['MINTmap_folder'] + 'OtherAnnotations.MINTmap_v1.txt'	
	MINTmap_MINTplates = config['EXECUTABLES']['MINTmap_folder'] + 'MINTplates/'
	
	results = []
	
	for jread in reads:
		sample_search = re.search('.*\/([a-zA-Z]{2,3})\_(\d{1,2})\_trimmed.*', jread)	
		if sample_search:
			outdir = sample_search.group(1) + "_" + sample_search.group(2)
			sample_folder =  folder + '/' + outdir + '/'
			if (os.path.isdir(sample_folder)):
				print ('\tMINTmap analysis for sample %s already exists' %outdir) 		
			else:
				#MINTmap.pl -f trimmedfastqfile [-p outputprefix] [-l lookuptable] [-s tRNAsequences] [-o tRFtypes] [-d customRPM] [-a assembly] [-j MINTplatesPath] [-h]
				fol = create_subfolder(outdir, folder)
				cmd = 'perl '+ MINTmap + ' -f %s -p %s -l %s -s %s -o %s -j %s' %(jread, sample_folder + outdir, MINTmap_table, MINTmap_tRNAseq, MINTmap_tRF, MINTmap_MINTplates) 

				results.append(sample_folder)

				# ToDOs: print into file MINTmap command
				print ('The following cmd is being executed at the shell: \n', cmd)
				try:
					subprocess.check_output(cmd, shell = True)
					print ('\tChecking sample %s for tRFs (MINTmap)' %outdir)
				except subprocess.CalledProcessError as err:
					print (err.output)
					sys.exit()
	
	return results		
###############

###############
def tRFs_analysis(path, count, reads, time_partial):
	##############################
    ####### Step: sRNAbench ######
	##############################
	print ("\n+ Run MINTmap: ")
	name_MINTmap_folder = str(count) + '.tRFs_MINTmap'
	MINTmap_folder = create_subfolder(name_MINTmap_folder, path)
	results = MINTmap(reads, MINTmap_folder)
	
	## timestamp
	time_partial = timestamp(time_partial)
###############
 
#################################################
######				MAIN					#####
#################################################
if __name__ == "__main__":

	start_time_total = time.time()
  
  	## control if options provided or help
	if len(sys.argv) > 1:
		print ("")
	else:
		help_options()
		exit()    	
	##
    	 
	print ("\n######## Starting Process ########")
	now = datetime.now()
	date_time = now.strftime("%m/%d/%Y, %H:%M:%S")
	print (date_time, "\n")
  
    ########################
    ## get configuration  ##
	########################
	configuration_path_file = os.path.abspath(argv[1])
	config = configparser.ConfigParser()
	config.read(configuration_path_file)

	print ("------ Get configuration options ------")
	print ("Reading configuration file: ", configuration_path_file, "\n")
	
	## dump config file to stdout
	f = open(configuration_path_file, 'r')
	print ("\n***********************************************")
	config = f.read()
	f.close()
	print (config)
	print ("\n***********************************************")
	
	## configuration options:
	merge_option = config["VARIABLES"]["merge_samples"]
	analysis_option = config["VARIABLES"]["option"]
	threads_option = config["VARIABLES"]["thread"]
    
    ####################
    ## create folder  ##	
 	####################
	print ("\n------ Create project folder ------")
	project_name = config["GENERAL"]["project"]
	try:
		folder_path = os.path.abspath(project_name)
	except OSError:  
		path = os.getcwd()  
		print ("The current working directory is %s", path)
		print ("Folder for project name: ",  project_name, " will be created here")
		folder_path = path.append(project_name)
	    
    # define the access rights
	access_rights = 0o755
	try:
		os.mkdir(folder_path, access_rights)
	except OSError:  
	   	print ("Directory %s already exists" % folder_path)
	else:  
		print ("Successfully created the directory %s " % folder_path)
	
	print ("")

    ## checking file(s) exists
	paired_end = False
	print ("------ Checking files provided ------")
	file_R1 = config["GENERAL"]["fastq_R1"]
	file_R2 = []
	if config.has_option("GENERAL","fastq_R2"):
		file_R2 = config["GENERAL"]["fastq_R2"]
		paired_end = True

	try:
		os.path.isfile(file_R1)
	except FileNotFoundError:
		print ("ERROR: No file R1 provided\n")
		print ("Exit")
		exit()
	else:
		print ("+ Folder for fastq R1 is readable and accessible")
    
	if paired_end:
		try:
		    os.path.isfile(file_R2)
		except FileNotFoundError:
			print ("+ Folder for fastq R2 does not exists")
			print ("+ Using single end option")
		else:
			print ("+ Folder for fastq R2 is readable and accessible")
			print ("+ Using paired-end option\n")
	else:
		print ("+ Using single end option")

    ## checking prefix provided
	print ("------ Checking prefix provided ------")
	prefix_list = config["VARIABLES"]["prefix"]

	if prefix_list == 'all':
		print ("+ All samples will be retrieved")
	else:
		prefix_list2 = prefix_list.split(",")
		prefix_list = prefix_list2
		for samples in prefix_list:
			print ("+ Samples with prefix %s will be retrieved" % samples)
	
	print ("")

	################################################
	print ("------ Starting pipeline ------")
	################################################
  	
	#####################################
    ####### Step1: select_samples #######
	#####################################
	all_list_R1 = []
	all_list_R2 = []

	print ("\n+ Select samples: ")
	if prefix_list == 'all':
		sample_listR1 = select_samples("all", file_R1)
		all_list_R1.extend(sample_listR1)

		if paired_end:
			sample_listR2 = select_samples("all", file_R2)
			all_list_R2.extend(sample_listR2)
	else:

		for samples in prefix_list:
			print ("\t+",samples,"samples")
			sample_listR1 = select_samples(samples, file_R1)
			all_list_R1.extend(sample_listR1)
			if paired_end:
				sample_listR2 = select_samples(samples, file_R2)
				all_list_R2.extend(sample_listR2)

	print ("")

	##################################
    ####### Step2: get_samples #######
	##################################
	sample_folder = []
	samplesR1 = []
	samplesR2 = []
	sample_folder = create_subfolder("1.samples", path=folder_path)

	### merge_option
	if merge_option == "YES":
		print ("\n+ Merge samples: ")
		samplesR1 = one_file_per_sample(all_list_R1, file_R1, sample_folder, read="R1")
		if paired_end:
			samplesR2 = one_file_per_sample(all_list_R2, file_R2, sample_folder, read="R2")
	else:
		print ("\n+ Retrieve samples:")
		samplesR1 = get_symbolic_link(all_list_R1, file_R1, sample_folder)
		if paired_end:
			samplesR2 = get_symbolic_link(all_list_R2, file_R2, sample_folder)
		
	## timestamp
	start_time_partial = timestamp(start_time_total)
	
	###############################
    ####### Step3: cutadapt #######
	###############################
	print ("\n+ Trimming samples: ")
	cutadapt_folder = create_subfolder("2.cutadapt", path=folder_path)
	trimmed_R1_return, trimmed_R2_return = cutadapt(samplesR1, samplesR2, sample_folder, cutadapt_folder)
	## timestamp
	start_time_partial = timestamp(start_time_partial)
		
	###############################
	####### Step: fastqjoin ######
	###############################
	folder_id = 3
	joined_read = []
	if paired_end:
		print ("\n+ Joining samples: ")
		name = str(folder_id) + ".fastqjoin"
		folder_id = folder_id + 1
		fastqjoin_folder = create_subfolder(name, path=folder_path)
		joined_read = fastqjoin(trimmed_R1_return, trimmed_R2_return, fastqjoin_folder)
		## timestamp
		start_time_partial = timestamp(start_time_partial)
	else:
		joined_read = trimmed_R1_return
		
	######################################
	####### Step: Small RNA analysis #####
	######################################
	if analysis_option == "isomiR":
		isomiR_analysis(folder_path, folder_id, joined_read, start_time_partial)

	elif analysis_option == "tRFs":
		tRFs_analysis(folder_path, folder_id, joined_read, start_time_partial)
				
	elif analysis_option == "all":
		isomiR_analysis(folder_path, folder_id, joined_read, start_time_partial)
		folder_id = folder_id + 1
		tRFs_analysis(folder_path, folder_id, joined_read, start_time_partial)
	else:
		print ("**ERROR: No valid analysis provided: ", analysis_option)
		print ("Please provide: isomiR|tRFs|all")
		

	print ("\n*************** Finish *******************")
	start_time_partial = timestamp(start_time_total)
	
	
#############################

