#usr/bin/env python

## useful imports
import time
import io
import os
import re
import subprocess
import sys
from datetime import datetime
import concurrent.futures

###############   
def gettime (start_time):
	total_sec = time.time() - start_time
	m, s = divmod(int(total_sec), 60)
	h, m = divmod(m, 60)
	
	return h, m, s


###############   
def create_subfolder (name, path):
    ## create subfolder  ##	
	subfolder_path = path + "/" + name
	access_rights = 0o755
	
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
def create_folder (path):
    ## create subfolder  ##	
	access_rights = 0o755
	
    # define the access rights
	try:
		os.mkdir(path, access_rights)
	except OSError:  
	   	print ("\tDirectory %s already exists" %path)
	else:  
		print ("\tSuccessfully created the directory %s " %path)
	
	print ("")
	return path
###############  

###############	
def timestamp (start_time_partial):
	h,m,s = gettime(start_time_partial)
	print ('--------------------------------')
	print ('(Time spent: %i h %i min %i s)' %(int(h), int(m), int(s)))
	print ('--------------------------------')
	return time.time()

############### 
def get_symbolic_link (sample_list, path_to_samples, directory):
	for samplex in sample_list:
		sample_path = path_to_samples + '/' + samplex
		cmd = 'ln -s %s %s' %(sample_path, directory)
		system_call(cmd)

	files2return = os.listdir(directory)
	return files2return

###############
def system_call(cmd):	
	## call system
	## send command	
	try:
		subprocess.check_output(cmd, shell = True)
		return ('OK')
	except subprocess.CalledProcessError as err:
		print (err.output)
		return ('FAIL')

###############	
def extract(fileGiven):
	print ("")
	#xtract(fileGiven, all=True)
	
###############
def sender(list_cmd, num_threads):	
	# We can use a with statement to ensure threads are cleaned up promptly
	with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as executor:
		# Start the load operations and mark each future with its URL
		commandsSent = { executor.submit(command_sender, commands): commands for commands in list_cmd }	
		for cmd2 in concurrent.futures.as_completed(commandsSent):
			details = commandsSent[cmd2]
			try:
				data = cmd2.result()
			except Exception as exc:
				print ('***ERROR:')
				print (string2send)
				print('%r generated an exception: %s' % (details, exc))
###############

###############
def command_sender(string2send):
	#print (string2send)
	try:
		subprocess.check_output(string2send, shell = True)
	except subprocess.CalledProcessError as err:
		print ('')
		
###############


