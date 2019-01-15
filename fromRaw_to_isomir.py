#usr/bin/env python

"""
    Author v0.1: Antonio Luna de Haro
    Author v1: Jose F Sanchez-Herrero

    This scripts runs sRNAtoolbox for a selection of RNAseq samples
!!!!    Make sure to activate miRTop! (activate-mirtop-0.3.11a0)    !!!!
"""
## useful imports
import sys
from sys import argv
import time
from datetime import datetime
import subprocess
import os
import re
import configparser


#####################
#### functions ######
#####################

###############
def select_samples (samples_prefix, path_to_samples):
    #Get all files in the folder "path_to_samples"    
    files = os.listdir(path_to_samples)
    sample_list = []
    for fastq in files:
    	#print "Test: ",fastq
        if fastq.startswith(samples_prefix) and 'merged' not in fastq: 
            if fastq.endswith('.gz'):
                sample_list.append(fastq[:-3]) 
            elif fastq.endswith('fastq'):
                sample_list.append(fastq)
            else:
                print "** ERROR: ", fastq, 'is a file that is neither in fastq.gz or .fastq format, so it is not included'
    
    non_duplicate_samples = list(set(sample_list))
    number_samples = len(non_duplicate_samples)
    print "\t\t- ", number_samples," samples selected from ", path_to_samples
    return sorted(non_duplicate_samples)
###############


###############
def gunziping (non_duplicate_samples, reverse = False, path_to_samples = "/imppc/labs/lslab/share/isomiR_project/mirXPlore/FASTQ"):
    path_to_samples = "/imppc/labs/lslab/aluna/Ilumina_NGS/150120_VHIO_data_FIS_SERUM_miRNAs/fastqR"
    if reverse:
        path_to_samples = path_to_samples + '2'
    else:
        path_to_samples = path_to_samples + '1'
    
    final_sample_list = []     
    for sample in non_duplicate_samples:
        sample = path_to_samples + sample
        if not os.path.isfile(sample):
            sample = sample + '.gz'
            print 'Using gunzip to -->', sample
            subprocess.check_output(['gunzip','-k', sample])
            sample = sample[:-3]
        final_sample_list.append(sample)
    return sorted(final_sample_list)
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
			#print "samplename_search", samplename_search.group(1)
			#print "samplename_search", samplename_search.group(2)
			#print "samplename_search", samplename_search.group(3)

			if samplename_search:
				path = samplename_search.group(1)
				sample = samplename_search.group(2)
				comonpart = path + "_" + sample + "_"
				commonname = path + "_" + sample + "_" + read + ".fastq"
				bigfilepath = directory + "/" + commonname
				bigfile_list.append(commonname)
					
				#print comonpart
				#print commonname
				#print bigfilepath
				
				for sampley in final_sample_list:
					if comonpart in sampley:
						#print "sample: ", sampley
						subsamples.append(path_to_samples + "/" + sampley)
						grouped_subsamples.append(sampley)
					
				if not os.path.isfile(bigfilepath) or os.stat(bigfilepath).st_size == 0:
					partsofsample = ' '.join(sorted(subsamples))
					cmd = 'cat %s >> %s' %(partsofsample, bigfilepath)
					##print "Merge command: ", cmd ## DUMP in file merge_info.txt
					try:
						print '\t+ %s created' %commonname
						subprocess.check_output(cmd, shell = True)
					except subprocess.CalledProcessError as err:
						print err.output
						sys.exit()
				else:
					print '\t + Sample %s is already merged' % commonname

	print 'There are' , len(bigfile_list) , 'samples after merging for read' , read, '\n'

	return bigfile_list
	
###############    

###############   
def cutadapt (list_R1, list_R2, path, out_path):
    trimed_R1 = []
    trimed_R2 = [] 

    for file_R1 in list_R1:
    
        sampleR1_search = re.search('(.*\/)([a-zA-Z]{2,3})(\d{5})(\d{2})(.*)', file_R1)
        if sampleR1_search:
            o_param = out_path + sampleR1_search.group(2) + '_' + sampleR1_search.group(4) + '_1.fastq'
            p_param = out_path + sampleR1_search.group(2) + '_' + sampleR1_search.group(4) + '_2.fastq'
        
            for file_R2 in list_R2:
                
                sampleR2_search = re.search('(.*\/)([a-zA-Z]{2,3})(\d{5})(\d{2})(.*)', file_R2)
                if sampleR2_search:
                    
                    if sampleR1_search.group(2) == sampleR2_search.group(2) and sampleR1_search.group(4) == sampleR2_search.group(4):
                        trimed_R1.append(o_param)
                        trimed_R2.append(p_param)
                        if not (os.path.isfile(o_param) or os.path.isfile(p_param)):
                            cmd = 'cutadapt -a TGGAATTCTCGGGTGCCAAGG -A GATCGTCGGACTGTAGAACTCTGAAC -o %s -p %s %s %s' %(o_param, p_param, file_R1, file_R2)
                            print 'The following cmd is being executed at the shell: \n', cmd
                            try:
                                print 'Trimming for %s sample' %(sampleR1_search.group(2) + '_' + sampleR1_search.group(4))
                                subprocess.check_output(cmd, shell = True)
                                print 'Adapters trimmed for the sample %s' %(sampleR1_search.group(2) + '_' + sampleR1_search.group(4))
                            except subprocess.CalledProcessError as err:
                                print err.output
                                sys.exit()                    
                        else:
                            print 'Sample %s is already trimmed' %(sampleR1_search.group(2) + '_' + sampleR1_search.group(4))
    return trimed_R1, trimed_R2
###############

###############     
def fastqjoin (trimed_R1, trimed_R2, error_param = '0.1'):
    out_path = '/imppc/labs/lslab/aluna/Ilumina_NGS/150120_VHIO_data_FIS_SERUM_miRNAs/fastqjoin/'
    joined_reads = []
    for read1 in trimed_R1:
        sample1 = read1.rpartition('/')[2].split('.')[0][:-2]
        for read2 in trimed_R2:
            sample2 = read2.rpartition('/')[2].split('.')[0][:-2]
            if sample1 == sample2:
                joined_reads.append(out_path + sample1 + '_join.fastq')
                if not (os.path.isfile(out_path + sample1 + '_join.fastq')):
                    
                    cmd = '/soft/bio/ea-utils-expressionanalysis-git/bin/fastq-join -p %s %s %s -o %s' %(error_param, read1, read2, out_path + sample1 + '_%.fastq')
                    print 'The following cmd is being executed at the shell: \n', cmd
                    
                    
                    try:
                        print 'Sample %s is joining the paired-end reads' %sample1
                        subprocess.check_output(cmd, shell = True)
                        
                    except subprocess.CalledProcessError as err:
                        print err.output
                        sys.exit()
                else: 
                    print 'Sample %s is already joined' %sample1
    return joined_reads
###############       
    
###############       
def sRNAbench (joined_reads, paired=True):
    outpath = '/imppc/labs/lslab/aluna/opt/sRNAtoolboxDB/out/'
    results = []
    
        
    for jread in joined_reads:
        if paired:
            outdir = jread.rpartition('/')[2].split('.')[0]
            finalpath = outpath + 'paired/' + outdir + '/'
        else:
            outdir = jread.rpartition('/')[2].split('.')[0][:-2]
            finalpath = outpath + 'single/' + outdir + '/'
        results.append(finalpath)
        print finalpath
        if not (os.path.isdir(finalpath)):
            cmd = 'java -jar /imppc/labs/lslab/aluna/opt/sRNAtoolboxDB/exec/sRNAbench.jar dbPath=/imppc/labs/lslab/aluna/opt/sRNAtoolboxDB/ input=%s output=%s microRNA=hsa  isoMiR=true plotLibs=true graphics=true plotMiR=true bedGraphMode=true writeGenomeDist=true chromosomeLevel=true chrMappingByLength=true ' %(jread, finalpath)
            print 'The following cmd is being executed at the shell: \n', cmd
            
            try:
                #print 'Searching isomiRs for %s' %outdir[:-5](.*\/)((([a-zA-Z]{2,3})_(C_)?(\d{1}|\d{2}))_[ACGT]{6}_(.{4})_R\d)_(.*)
                subprocess.check_output(cmd, shell = True)
                
            except subprocess.CalledProcessError as err:
                print err.output
                sys.exit()
        else: 
            print 'Sample %s has already been serached for isomiRs' %outdir[:-5]
    return results
###############   
    
###############   
def miRTop (results):
    outpath = '/imppc/labs/lslab/aluna/q/'
    gtfs = []
    for folder in results:
        outdir = '/'.join(folder.split('/')[-3:-1])
        print outdir
        gtfs.append(outpath + outdir + '/')
        if not (os.path.isdir(outpath + outdir)):
            cmd = 'mirtop gff --sps hsa --hairpin /imppc/labs/lslab/aluna/opt/sRNAtoolboxDB/libs/hairpin.fa --gtf /imppc/labs/lslab/aluna/opt/hsa.gff3 --format srnabench -o %s  %s' %(outpath + outdir + '/', folder)

            print 'The following cmd is being executed at the shell: \n', cmd
            
            try:
                print 'Creating isomiRs gtf file for %s' %outdir
                subprocess.check_output(cmd, shell = True)
                
            except subprocess.CalledProcessError as err:
                print err.output
                sys.exit()
        else: 
            print 'Sample %s has already an isomiRs gtf file' %outdir
    return gtfs
    #return results
###############   
    
###############   
def miRTop_stats (gtfs):
    for gtf in gtfs:
        gtfname = gtf.split('/')[-2]
        gtfile = gtf + gtfname + '.gtf'
        outdir = gtf + 'Stats'
        if not (os.path.isdir(outdir)):
            cmd = 'mirtop stats -o %s %s' %(outdir, gtfile) 
            print 'The following cmd is being executed at the shell: \n', cmd
            try:
                print 'Creating isomiRs stats for %s' %gtfname
                subprocess.check_output(cmd, shell = True)
            except subprocess.CalledProcessError as err:
                print err.output
                sys.exit()
        else: 
            print 'Sample %s has already isomiRs stats' %gtfname
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
	   	print ("Directory %s already exists" % subfolder_path)
	else:  
		print ("Successfully created the directory %s " % subfolder_path)
	
	print ""
	return subfolder_path
###############   
    
###############   
def timestamp (start_time_partial):
	h,m,s = getime(start_time_partial)
	print '---------------------------------------------------'
	print '(Time spent: %i h %i min %i s)' %(int(h), int(m), int(s))
	print '---------------------------------------------------'
	return time.time()
###############   
    
def help_options():
	print "\n#################################################"
	print "Usage:\npython", os.path.abspath(argv[0]),"config_file.txt "
	print ""
	print "This script... [Write a description]"
	print ""
	print "Configuration file includes general and detailed information for the project:"
	print "For further details check: (...)"
	print "--"
	print "[GENERAL]"
	print "fastq_R1 = /path/to/file/fastqR1"
	print "fastq_R2 = /path/to/file/fastqR2"
	print "project = project_name"
	print "prefix = prefix_name"
	print ""
	print "[EXECUTABLES]"
	print "cutadapt = /path/to/cutadapt/bin"
	print "--"
	print "#################################################\n\n"
	

#################################################
######				MAIN					#####
#################################################
if __name__ == "__main__":
	start_time_total = time.time()
  
  	## control if options provided or help

	if len(sys.argv) > 1:
		print ""
	else:
		help_options()
		exit()
    	
    	 
	print "\n######## Starting Process ########"
	now = datetime.now()
	date_time = now.strftime("%m/%d/%Y, %H:%M:%S")
	print date_time, "\n"
  
    ########################
    ## get configuration  ##
	########################
	configuration_path_file = os.path.abspath(argv[1])
	config = configparser.ConfigParser()
	config.read(configuration_path_file)

	print "------ Get configuration options ------"
	print "Reading configuration file: ", configuration_path_file, "\n"
    
    ####################
    ## create folder  ##	
 	####################
	print "------ Create project folder ------"
	project_name = config["GENERAL"]["project"]
	try:
		folder_path = os.path.abspath(project_name)
	except OSError:  
		path = os.getcwd()  
		print "The current working directory is %s", path
		print "Folder for project name: ",  project_name, " will be created here"
		folder_path = path.append(project_name)
	    
    # define the access rights
	access_rights = 00755
	try:
		os.mkdir(folder_path, access_rights)
	except OSError:  
	   	print ("Directory %s already exists" % folder_path)
	else:  
		print ("Successfully created the directory %s " % folder_path)
	
	print ""
	
    #####################################################
    ## 	Reading arguments from the configuration file  ##
    #####################################################

    ## checking file(s) exists
	print "------ Checking files provided ------"
	file_R1 = config["GENERAL"]["fastq_R1"]
	file_R2 = config["GENERAL"]["fastq_R2"]

	try:
		os.path.isfile(file_R1)
	except FileNotFoundError:
		print "ERROR: No file R1 provided\n"
		print "Exit"
		exit()
	else:
		print "+ Folder for fastq R1 is readable and accessible"
    
	try:
	    os.path.isfile(file_R2)
	except FileNotFoundError:
		print "+ Folder for fastq R2 does not exists"
		print "+ Using single end option"
	else:
		print "+ Folder for fastq R2 is readable and accessible"
		print "+ Using paired-end option\n"
		paired_end = True

    ## checking prefix provided
  	print "------ Checking prefix provided ------"
	prefix_list = config["GENERAL"]["prefix"]

	if prefix_list == 'all':
		print "+ All samples will be retrieved"
	else:
		prefix_list2 = prefix_list.split(",")
		prefix_list = prefix_list2
		for samples in prefix_list:
			print ("+ Samples with prefix %s will be retrieved" % samples)
	
	print ""
		

	################################################
  	print "------ Starting pipeline ------"
	################################################
  	
	#####################################
    ####### Step1: select_samples #######
	#####################################
	all_list_R1 = []
	all_list_R2 = []

	print "+ Select samples: "
	for samples in prefix_list:
		print "\t+",samples,"samples"
		sample_listR1 = select_samples(samples, file_R1)
		all_list_R1.extend(sample_listR1)
		if paired_end:
			sample_listR2 = select_samples(samples, file_R2)
			all_list_R2.extend(sample_listR2)

	####################################
    ####### Step2: merge_samples #######
	####################################
	print "+ Merge samples: "
	merge_folder = create_subfolder("1.merge", path=folder_path)
	mergeR2 = []
	mergedR1 = one_file_per_sample(all_list_R1, file_R1, merge_folder, read="R1")
	if paired_end:
		mergedR2 = one_file_per_sample(all_list_R2, file_R2, merge_folder, read="R2")
	## timestamp
	start_time_partial = timestamp(start_time_total)
	
	exit()

	###############################
    ####### Step3: cutadapt #######
	###############################
	cutadapt_folder = create_subfolder("2.cutadapt", path=folder_path)
	trimed_R1, trimed_R2 = cutadapt(mergedR1, mergedR2, merge_folder, cutadapt_folder)

	## timestamp
	start_time_partial = timestamp(start_time_partial)
	
	###############################
    ####### Step4: fastqjoin ######
	###############################
	fastqjoin_folder = create_subfolder("3.fastqjoin", path=folder_path)
	joined_reads =fastqjoin(trimed_R1, trimed_R2)
	## timestamp
	start_time_partial = timestamp(start_time_partial)

	###############################
    ####### Step5: sRNAbench ######
	###############################
	sRNAbench_folder = create_subfolder("4.sRNAbench", path=folder_path)
	if type_reads == 'paired':    
		results = sRNAbench(joined_reads)
	else:
		results = sRNAbench(trimed_R1, paired=False)
	## timestamp
	start_time_partial = timestamp(start_time_partial)

	#############################
    ####### Step6: miRTop #######
	#############################
	miRTop_folder = create_subfolder("5.miRTop", path=folder_path)
	gtfs=miRTop(results)
	h,m,s = getime(start_time_partial)
	## timestamp
	start_time_partial = timestamp(start_time_partial)

	###################################
    ####### Step6: miRTop_stats #######
	###################################
	miRTop_stats(gtfs)
	print "Finish"
	start_time_partial = timestamp(start_time_total)

