#usr/bin/env python

## useful imports
import time
import io
import os
import re
import subprocess
import sys
from sys import argv
from datetime import datetime
from io import open


## https://www.biostars.org/p/310096/

#####################
#### functions ######
#####################

def help_options():
	print ("\n#######################################################################")
	print ("  NAME: RNA_biotypes")
	print ("  VERSION: 0.1")
	print ("  AUTHORS: Jose F Sanchez-Herrero (v1).")
	print ("           Copyright (C) 2018-2019 Lauro Sumoy Lab, IGTP, Spain")
	print ("#########################################################################")
	print ("\nDESCRIPTION:")
	print ("- This script is a pipeline for the identification of RNA biotypes.\n")
	print ("USAGE:\npython3", os.path.abspath( argv[0] ),"STAR_executable genomeDir cpus outName fastq_R1 [fastq_R2]")
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
def gettime (start_time):
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
	h,m,s = gettime(start_time_partial)
	print ('--------------------------------')
	print ('(Time spent: %i h %i min %i s)' %(int(h), int(m), int(s)))
	print ('--------------------------------')
	return time.time()
############### 


#################################################
######				MAIN					#####
#################################################
if __name__ == "__main__":

	start_time_total = time.time()
  
  	## control if options provided or help
	if len(sys.argv) > 5:
		print ("")
	else:	
		help_options()
		exit()    	

	## Get arguments provided
	STAR_exe = argv[1]
	genomeDir = argv[2]
	cpus =  argv[3]
	outName = argv[4]
	read_file1 = argv[5]
	read_file2 = []
	
	if (argv[5]):
		read_file2 = argv[6]

    	 
	print ("\n######## Starting Process ########")
	now = datetime.now()
	date_time = now.strftime("%m/%d/%Y, %H:%M:%S")
	print (date_time, "\n")
	print ("+ Map and obtain mapping reads using STAR:")

	cmd = "%s --genomeDir %s --runThreadN %s --readFilesIn %s" %(STAR_exe, genomeDir, cpus, read_file1)
	if read_file2:
		cmd = cmd + " %s" %read_file2
	
	## all this options and parameters have been obtained from https://www.encodeproject.org/rna-seq/small-rnas/
	cmd = cmd + " --outFilterMultimapNmax 20 --alignIntronMax 1 --outFilterMismatchNoverLmax 0.03 --outFilterScoreMinOverLread 0 "
	cmd = cmd + "--outFilterMatchNminOverLread 0 --outFilterMatchNmin 16 --outSAMheaderHD @HD VN:1.4 SO:coordinate --outSAMunmapped Within "
	cmd = cmd + "--outSAMtype BAM SortedByCoordinate --genomeLoad NoSharedMemory --quantMode GeneCounts --alignSJDBoverhangMin 1000 "
	cmd = cmd + "--outFileNamePrefix %s" %outName
	
	print ("+ Command:")
	print ("\n#############################")
	print (cmd)
	print ("#############################\n")

	try:
		subprocess.check_output(cmd, shell = True)
	
	except subprocess.CalledProcessError as err:
	
		print (err.output)
		sys.exit()

