#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez										##
## Copyright (C) 2019-2020 Lauro Sumoy Lab, IGTP, Spain	##
##########################################################
"""
Provides configuration for the pipeline.
"""
## useful imports
import os
import io
import sys
import re
import shutil
from io import open
from sys import argv
import subprocess
from termcolor import colored
from distutils.version import LooseVersion
import pkg_resources

## import my modules
from HCGB import functions
from XICRA.config import extern_progs

################
## Software
################

##################
def get_exe(prog, Debug=False, Return_Version=False):
	"""Return absolute path of the executable program requested.

	Given a program name it returns its executable to be called. It has to fulfilled a minimum version specified.

	:param prog: Software name
	:type prog: string
	:returns: Absolute path for the executable requested
	:warning: if no executable available in system ``$PATH`` or not fulfilling the expected version.

	.. attention:: Be aware of Copyright

		The code implemented here was retrieved and modified from ARIBA (https://github.com/sanger-pathogens/ariba)

		Give them credit accordingly.

	"""
	exe = ""
	if prog in os.environ: 
		exe = os.environ[env_var] ## python environent variables
	else:
		exe = extern_progs.return_defatult_soft(prog) ## install in the system

	## get paths
	exe_path_tmp = my_which(exe)

	## debug message
	if (Debug):
		print(colored("** Debug: exe: %s" %exe,'yellow'))
		print(colored("** Debug: exe_path_tmp: %s" %exe_path_tmp,'yellow'))

	## get min_version
	min_version = extern_progs.return_min_version_soft(prog)
	
	## debug message
	if (Debug):
		print(colored("** Debug: min_version: %s" %min_version,'yellow'))

	## return if not available
	## no min version available
	if min_version == 'na':
		if exe_path_tmp:
			if (Return_Version):
				return (exe_path_tmp[0], '') ## return first item
			else:
				return (exe_path_tmp[0]) ## return first item
	
	## not installed in path
	if exe_path_tmp is None or len(exe_path_tmp) == 0:
		if (Return_Version):
			print(colored("\n**ERROR: Software %s could not be found." % prog,'red'))
			exit()
			return('ERROR', 'n.a.')
		else:
			print(colored("\n**ERROR: Software %s could not be found." % prog,'red'))
			exit()
			return('ERROR')
		

	## Loop for all possibilities
	for p in exe_path_tmp:
		prog_ver = get_version(prog, p, Debug=Debug)

		if (Debug):
			print (colored("** Debug: Software: %s\nPath: %s\nVersion: %s" %(prog, p, prog_ver), 'yellow'))

		if (prog_ver == 'n.a.'):
			continue

		if LooseVersion(prog_ver) >= LooseVersion(min_version):
			if (Return_Version):
				return (p, prog_ver)
			else:
				return (p)


	print(colored("\n**ERROR: Software %s version smaller than minimum version expected %s." %(prog,min_version),'red'))
	exit()

################
def access_check(fn, mode=os.F_OK | os.X_OK):
	"""Check executable permission

	This function checks whether a given path is a folder or file and if it is 
	executable and accessible. It also works if a java jar file provided.

	:param fn: Absolute path file
	:param mode: Value to pass as the mode parameter of access()

	:type fn: string
	:type mode: string

	`mode` defaults to:

		- os.F_OK: Value to pass as the mode parameter of access() to test the existence of path.

		- os.X_OK: Value to include in the mode parameter of access() to determine if path can be executed.

	.. attention:: Be aware of Copyright

		The code implemented here was retrieved and modified from shutil (https://github.com/python/cpython/blob/master/Lib/shutil.py).

		Give them credit accordingly.

		We modified the code to work if java jar files provided.
	"""
	## the original code belongs to shutil, slightly modified here
	# https://github.com/python/cpython/blob/master/Lib/shutil.py

	#if os.path.isdir(fn):
	#	return False

	if os.path.exists(fn):
		if fn.endswith('.jar'):
			return True

		if os.access(fn, mode):
			return True

#################
def my_which(cmd):
	"""Return the absolute path to the executable

	Given a command return the absolute path(s), if any.

	:param cmd: Software command name 
	:returns: List of absolute paths(s) of the given command.

	.. attention:: Be aware of Copyright

		The code implemented here was retrieved and modified from shutil (https://github.com/python/cpython/blob/master/Lib/shutil.py).

		Give them credit accordingly.

		We modified the code to return multiple paths in a list if available different installed binaries in $PATH.

	"""
	# If we're given a path with a directory part, look it up directly rather
	# than referring to PATH directories. This includes checking relative to the
	# current directory, e.g. ./script
	if os.path.dirname(cmd):
		if access_check(cmd):
			return cmd
		return None

	use_bytes = isinstance(cmd, bytes)

	path=None

	if path is None:
		path = os.environ.get("PATH", None)

	if path is None:
		try:
			path = os.confstr("CS_PATH")
		except (AttributeError, ValueError):
			# os.confstr() or CS_PATH is not available
			path = os.defpath
		# bpo-35755: Don't use os.defpath if the PATH environment variable is
		# set to an empty string

	# PATH='' doesn't match, whereas PATH=':' looks in the current directory
	if not path:
		return None

	if use_bytes:
		path = os.fsencode(path)
		path = path.split(os.fsencode(os.pathsep))
	else:
		path = os.fsdecode(path)
		path = path.split(os.pathsep)

	# On other platforms you don't have things like PATHEXT to tell you
	# what file suffixes are executable, so just pass on cmd as-is.
	files = [cmd]

	return_paths = [] ## modification
	seen = set()
	for dir in path:
		normdir = os.path.normcase(dir)
		#print ("Normdir: ", normdir)
		if not normdir in seen:
			seen.add(normdir)
			for thefile in files:
				name = os.path.join(dir, thefile)
				#print ("Name: ", name)
				if access_check(name):
					## return (name) ## previously, it would only return the first item
					return_paths.append(name) ## modification

	if (len(return_paths) >= 1):
		return return_paths
	else:
		return None

##################
def get_version(prog, path, Debug=False):
	"""Get version of software

	Given a program name and expected path, tries to determine its version.

	:param prog: Program name
	:param path: Absolute path
	:param Debug: True/False

	:type prog: string
	:type path: string 
	:type Debug: bool

	:returns: String containing version. Returns NA message if no found and raises attention error message.

	.. attention:: Be aware of Copyright

		The code implemented here was retrieved and modified from ARIBA (https://github.com/sanger-pathogens/ariba)

		Give them credit accordingly.
	"""

	## read dependencies information
	dependencies_pd = extern_progs.read_dependencies()

	## get information for prog
	regex = re.compile(dependencies_pd.loc[prog, 'get_version'])
	args = dependencies_pd.loc[prog, 'version_cmd']
	cmd = path + ' ' + args

	## debug messages
	if (Debug):
		print(colored("** Debug: regex: %s" %regex,'yellow'))
		print(colored("** Debug: args: %s" %args, 'yellow'))

	if prog == 'sRNAbench' or prog == "miraligner":
		java_bin = get_exe('java', Debug=Debug)
		java_jar = java_bin + ' -jar ' + path + ' ' + args
		cmd_output = subprocess.Popen(java_jar, shell=True, stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()
	else:
		cmd_output = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()

	## decode command
	cmd_output = functions.main_functions.decode(cmd_output[0]).split('\n')[:-1] + functions.main_functions.decode(cmd_output[1]).split('\n')[:-1]

	## debug messages
	if (Debug):
		print(colored("** Debug: cmd_output:\n %s" %cmd_output,'yellow'))

	## retrieve version information
	for line in cmd_output:
		hits = regex.search(line)
		if (Debug):
			print (hits)
		if hits:
			return hits.group(1)

	if Debug:
		print (colored('Attention: I tried to get the version of ' + prog + ' with: "' + cmd + '" and the output didn\'t match this regular expression: "' + regex.pattern + '"', 'red'))

	return("n.a.")

def check_dependencies(Debug):
	"""
	Check if available the different software required for ``XICRA`` execution.
	
	Using the function :func:`XICRA.config.extern_progs.read_dependencies` the 
	information for all the dependencies is retrieved from file :file:`XICRA/config/software/dependencies.csv`.
	
	For each software, the function :func:`XICRA.config.set_config.get_exe` retrieves
	whether it is installed in the system or not and its version and it is check using 
	:func:`XICRA.config.set_config.check_install`. 
	
	:param Debug: True/False for debugging messages 
	:type Debug: boolean
	
	:returns: Print messages and information	
	"""
	
	## read dependencies information
	dependencies_pd = extern_progs.read_dependencies()

	for soft, row in dependencies_pd.iterrows():
		soft_name = row['soft_name']
		min_version = row['min_version']
		
		## debug messages
		if (Debug):
			print ("Software:", soft)
			print ("Soft name:", soft_name)
			print ("Min_Version:", min_version)
		
		## get executable
		(soft_path, installed) = get_exe(soft, Debug=Debug, Return_Version=True)
		
		## debug messages
		if (Debug):
			print ("Soft Path: ", soft_path)
			print ("Version installed:", installed)
			
		## check if installed
		if (min_version == 'na'):
			message = check_install_module('ok', soft_name, min_version, 'Software')
		else:
			message = check_install_module(installed, soft_name, min_version, 'Software')

		if (message == 'OK'):
			continue
		else:
			print ("+ Please install manually software: ", soft_name, " to continue with XICRA\n\n")
	
################
## Python
################
def get_python_packages(Debug):
	"""
	Retrieves the version of the python packages installed in the system.

	It retrieves the dependencies name conversion from file :file:`XICRA/config/python/module_dependencies.csv`
	using function :func:`XICRA.config.extern_progs.file_list` and :func:`XICRA.scripts.functions.main_functions.get_data`.
	For each module it retrieves the package version installed in the system using 
	:func:`XICRA.config.set_config.check_package_version`.	

	:returns: Dictionary containing for each python module (key) the installed version (value).
	"""

	## get import names for packages:
	## some modules do not have the same name when install from pip and called from import
	file_module_dependecies = extern_progs.file_list("python_requirements")
	module_dependencies = functions.main_functions.file2dictionary(file_module_dependecies, ',')

	my_packages_installed = {}
	for each in module_dependencies:
		module_name = module_dependencies[each]
		installed = check_package_version(each, Debug) ## check version installed in system
		my_packages_installed[each] = installed

	return (my_packages_installed)

##################
def check_python_packages(Debug):
	"""
	This functions checks whether the packages installed in the system fulfilled the 
	minimum version specified in the configuration folder. 

	It uses function :func:`XICRA.config.set_config.get_python packages` to
	retrieve the version of the python packages installed in the system. Then it uses
	:func:`XICRA.config.extern_progs.min_python_module_version` to retrieve the minimum
	version specified. It compares them using function :func:`XICRA.config.set_config.check_install_module`.

	:param Debug: True/False for debugging messages
	:type Debug: boolean

	:returns: Print messages if packages are installed.
	"""
	## get python packages installed
	my_packages_installed = get_python_packages(Debug)

	## debug messages
	if (Debug):
		print ("my_packages_installed :: ")
		print (my_packages_installed)

	## min versions for packages
	my_packages_requirements = extern_progs.min_python_module_version()

	## debug messages
	if (Debug):
		print ("my_packages_requirements")
		print (my_packages_requirements)

	## check each package
	for each in my_packages_requirements:
		## get min version
		min_version = my_packages_requirements[each]

		## get version installed in system
		installed = my_packages_installed[each] 

		## debug messages
		if (Debug):
			print ("Module:", each)
			print ("Min_Version:", min_version)
			print ("Version installed:", installed)
			
		## check if installed
		message = check_install_module(installed, each, min_version, 'Module')

		if (message == 'OK'):
			continue
		else:
			print ("+ Please install manually package: ", each, " to continue with XICRA\n\n")
			#print ("pip install %s" %each)

################
def check_package_version(package, Debug):
	"""
	Retrieve python package version installed

	This is a modification of the original code from ARIBA (https://github.com/sanger-pathogens/ariba). 
	It basically uses pkg_resources.get_distribution(), pkg_resources.resource_filename() or imports module
	and retrieves version from __version__ variable.

	:param package: Python package name 
	:param Debug: True/False for debugging messages

	:type package: string
	:type Debug: boolean

	:returns: Version retrieved

	.. attention:: Be aware of Copyright

		The code implemented here was retrieved and modified from ARIBA (https://github.com/sanger-pathogens/ariba)

		Give them credit accordingly.
	"""

	try:
		version = pkg_resources.get_distribution(package).version
		if (Debug):
			print ("Method: pkg_resources.get_distribution(package).version")
	except:
		try:
			exec('import ' + package)
			version = eval(package + '.__version__')
			if (Debug):
				print ("Method: exec('import ' + package); version = eval(package + '.__version__')")
		except:
			try:
				if (Debug):
					print ("Method: pkg_resources.resource_filename(package, 'version.py')")
				version = pkg_resources.resource_filename(package, 'version.py')
			except:
				version = 'n.a.'
	if (Debug):
		print ('Package:', package)
		print ('Version:', version)
	return(version)
################

################
## R
################
def get_R_packages():
	dep_file = os.path.abspath(os.path.join(os.path.dirname( __file__ ), 'R', 'R_dependencies.csv'))
	dep_file_data = functions.main_functions.get_data(dep_file, ',', 'index_col=0')
	return (dep_file_data)

################
def check_R_packages(Debug):
	
	packages = get_R_packages()
	check_install_system = os.path.abspath(os.path.join(os.path.dirname( __file__ ), 'R', 'check_install_system.R'))
	R_script_exe = get_exe('Rscript')
	
	for index,row in packages.iterrows():
		## debugging messages
		if Debug:
			print ('\n+ Check package: ', index)
			print('+ Source: ', row['source'])
		
		cmd_check = R_script_exe + ' ' + check_install_system + ' -l ' + index
		code = functions.system_call_functions.system_call(cmd_check, message=False, returned=False)
		if (code=='OK'):
			check_install_module('1', index, '0', 'package')
		else:
			check_install_module('0', index, '1', 'System package')

			## check if installed in path
			cmd_check_path = R_script_exe + ' ' + check_install_path + ' -l ' + index + ' -p ' + install_path
			code2 = functions.system_call_functions.system_call(cmd_check_path, message=False, returned=False)

			if (code2):
				check_install_module('1', index, '0', 'Install path package')
			else:
				check_install_module('0', index, '1', 'Install path package')
				print ("Please install module %s manually to continue with XICRA" %index)
				
		
################
## Miscellaneous
################
def print_module_comparison(module_name, message, color, tag):
	"""
	Creates print message for a given module, version, message

	:param module_name: Name of the module
	:param message: Message to include in the print message: OK | FAILED | NOT FOUND
	:param color: Print message color: green | orange | red
	:param tag: Tag to include: Module, package, software

	:type module_name: string
	:type message: string 
	:type color: string
	:type tag: string

	:returns: Print message
	"""
	print (colored("{:.<15}{:.>15}".format("%s: %s" %(tag, module_name), "[ %s ]" %message), color))

#########

##################

def check_install_module(installed, module_name, min_version, tag):
	"""
	Checks modules installation 

	Checks whether a module is installed and fulfilling requirements. 	
	It prints messages using :func:`XICRA.config.set_config.print_module_comparison`.

	:param installed: Version string of the module installed.
	:param module_name: Module name
	:param min_version: Version string for the minimum version required

	:type installed: string
	:type module_name: string 
	:type min_version: string
	"""
	 ## Not installed
	if (installed == 'n.a.' or not installed):
		message = 'NOT FOUND'
		color = 'red'
		print_module_comparison(module_name, message, color, tag)

	# check version
	elif LooseVersion(installed) >= LooseVersion(min_version):
		message = 'OK'
		color = 'green'
		print_module_comparison(module_name, message, color, tag)

	# check version
	elif (installed == 'ok'):
		message = 'OK'
		color = 'green'
		print_module_comparison(module_name, message, color, tag)

	else:
		message = 'FAILED'
		color = 'yellow'
		print_module_comparison(module_name, message, color, tag)

	## return message
	return (message)

#########	
