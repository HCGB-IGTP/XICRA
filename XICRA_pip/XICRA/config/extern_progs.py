#!/usr/bin/env python3
##########################################################
## this modules is an idea from ARIBA					##
## (https://github.com/sanger-pathogens/ariba)			##
## give credit to them appropiately						##
##														##
## Jose F. Sanchez										##
## Copyright (C) 2019-2020 Lauro Sumoy Lab, IGTP, Spain	##
##########################################################
"""
Provides external programs details and configuration.
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
import pandas as pd
from termcolor import colored
from distutils.version import LooseVersion
import pkg_resources

## import my modules
from HCGB import functions
from XICRA.config import set_config

####################################################################
def file_list(wanted_data):
	"""
	Retrieves information of additional files under folder ``XICRA/config``.

	Using :func:`HCGB.functions.main_functions.get_fullpath_list` retrieves absolute
	path for file of interest.

	:param wanted_data: name for file
	:type wanted_data: string

	:returns: Absolute path for file wanted

	"""
	config_folder = os.path.dirname(os.path.realpath(__file__))
	listOffiles = functions.main_functions.get_fullpath_list(config_folder, False)
	
	for f in listOffiles:
		name = os.path.splitext(os.path.basename(f))[0]
		if (name == wanted_data):
			return (f)


##################
## Software
##################
def read_dependencies():
	"""Returns a dictionary containing the executable name for each software.

	It uses :func:`XICRA.config.extern_progs.file_list` to retrieve absolute path
	for file :file:`XICRA/config/software/dependencies.csv`. It then reads csv into pandas
	dataframe using :func:`XICRA.scripts.functions.main_functions.get_data` and returns it.	
	"""

	## read from file: prog2default.csv
	dependencies_file = file_list("dependencies")
	return(functions.main_functions.get_data(dependencies_file, ',', 'index_col=0'))

#######################
def return_defatult_soft(soft):
	"""Returns default name for a given software name

	For some software we provide a shorter name for the software. Here, we read file
	:file:`XICRA/config/software/dependencies.csv` using :func:`XICRA.config.extern_progs.read_dependencies`
	and retrieve original software name.

	:param soft: Software name
	:type soft: string

	:returns: String with original software name

	.. seealso:: This function depends on other XICRA functions:

		- :func:`XICRA.config.extern_progs.read_dependencies`

	"""
	dependencies_df = read_dependencies()
	return(dependencies_df.loc[soft,"soft_name"])

##################
def return_min_version_soft(soft):
	"""Retrieve version for a given software

	Retrieves minimum version for the software of interest stored in :file:`XICRA/config/software/dependencies.csv`.
	It reads file using :func:`XICRA.config.extern_progs.read_dependencies`
	and retrieve minimum version required.

	:param soft: Software name
	:type soft: string

	:returns: String with minimum version

	.. seealso:: This function depends on other XICRA functions:

		- :func:`XICRA.config.extern_progs.read_dependencies`	
	"""
	dependencies_df = read_dependencies()
	return(dependencies_df.loc[soft,"min_version"])
##################

##################
def print_dependencies():
	"""

	"""
	progs = {}
	depencencies_pd = read_dependencies()
	for prog in depencencies_pd:
		#print (prog)
		prog_exe = set_config.get_exe(prog)
		#print (prog + '\t' + prog_exe)
		prog_ver = get_version(prog, prog_exe)
		progs[prog] = [prog_exe, prog_ver]

	df_programs = pd.DataFrame.from_dict(progs, orient='index', columns=('Executable path', 'Version'))
	df_programs = df_programs.stack().str.lstrip().unstack()
	pd.set_option('display.max_colwidth', -1)
	pd.set_option('display.max_columns', None)
	print (df_programs)

##################
### Python packages
##################
def min_python_module_version():
	"""Returns a dictionary containing minimum version for each python package.

	Reads information from :file:`XICRA/config/python/python_requirements.csv`.

	:returns: dictionary
	"""
	## read from file: prog2default.csv
	python_modules = file_list("python_requirements")
	package_min_versions = functions.main_functions.file2dictionary(python_modules, ",")

	return(package_min_versions)
##################

##################
def return_min_version_python_package(package):
	"""
	Retrieves minimum version requirement for the given package.

	It retrieves the requirements using :func:`XICRA.config.extern_progs.min_python_module_version`
	and returns the given package requested minimun version.

	:param package:  
	:type package: string	
	:returns: Minimum version package (string)

	"""
	version_package = min_python_module_version()
	return (version_package[package])

##################
def print_package_version():
	"""
	Prints the package version required by ``XICRA``

	It retrieves the requirements using :func:`XICRA.config.extern_progs.min_python_module_version`
	and prints them using function :func:`XICRA.config.set_config.print_module_comparison`.

	:returns: Print messages
	"""
	my_packages = min_python_module_version()
	for each in my_packages:
		set_config.print_module_comparison(each, my_packages[each], 'green')
