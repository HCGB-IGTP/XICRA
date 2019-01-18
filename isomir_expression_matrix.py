#usr/bin/env python

"""
	AUTHOR:
    Antonio Luna de Haro (v0.1) & Jose F Sanchez-Herrero (v1)
	Copyright (C) 2018-2019 Lauro Sumoy Lab, IGTP, Spain

    DESCRIPTION:
    This scripts generates expression matrices for miRNA samples
    using the output generated from miRTop.

    
	LICENSE:
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
## useful imports
import sys
from sys import argv
import os

#####################
#### functions ######
#####################

def help_options():
	print "\n#################################################"
	print "isomir_expression_matrix\n"
	print "Version 1.0.0"
	print "Copyright (C) 2018-2019 Lauro Sumoy Lab, IGTP, Spain"
	print ""
	print "Usage:\npython", os.path.abspath(argv[0]),"folder4samples outfolder"
	print ""
	print "\nDescription:"
	print "This scripts generates expression matrices for miRNA samples using the output generated from miRTop."
	print ""
	print "\nParameters:"
	print "folder4samples: folder containing a folder for each sample of interest. Each subfolder should contian miRTop.gtf file"
	print "\tex: /path/to/folder_provided/sample_1/mirtop.gtf"
	print "outfolder: name for output folder\n"
	print "Citation:"
	print "[to add citation]"
	print ""
	print "See [ http://website ] for full documentation"
	print ""
	print ""
	print "#################################################\n\n"
###############

###############   
def getfiles_isomir(samples_prefix):
    files_s = os.listdir(samples_prefix)
    sample_s_list = []
    
    for singles in files_s:
        if singles.startswith(samples_prefix): 
           sample_s_list.append(files_s + "/mirtop.gtf" )
    
    isomirs_s = list( set(sample_s_list) )
    return sorted(isomirs)
    
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
	
	print ""
	return subfolder_path
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
    

if __name__ == "__main__":

	## argument provided: folder containing subfolders for each 
	## samples that contain miRTop output stats
	
	## control if options provided or help
	if len(sys.argv) > 1:
		print ""
	else:
		help_options()
		exit()    	

	##	
	prefix1 = argv[1]
	outfolder = argv[2]
	
	####################
    ## create folder  ##	
 	####################
	print "------ Create project folder ------"
	try:
		folder_path = os.path.abspath(outfolder)
	except OSError:  
		path = os.getcwd()  
		print "The current working directory is %s", path
		print "Folder for project name: ",  outfolder, " will be created here"
		folder_path = path.append(outfolder)
	else:
		print "Folder %s already exists" %outfolder
	print ""
	
	############################
	## retrieve miRTop files
	print "------ Get miRTop files generated ------"
	file_list = getfiles_isomir(prefix1)
	print ""

	##########################
	## Generate expression matrix
	filenames = []
	print "------ Generate expression matrix ------"
	for gtffile in file_list:

		sample = gtffile.rpartition('/')[-1][:-4]
		print "\t Parsing sample ", sample
		
		## parse gtf file
		sample_dict = parse_gtf(gtffile)

		## get filename
		filename = outfolder + '/' + sample
		filenames.append(filename)

		## create matrix
		make_table(sample_dict, filename)
		print '\tExpression matrix created for ', sample   

    
