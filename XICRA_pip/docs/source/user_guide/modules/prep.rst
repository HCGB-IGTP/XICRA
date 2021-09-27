.. ############################
.. _prep-description:

prep
====

This module prepares fastq files for later usage.

It initially checks the length of the name and advises the user to rename samples 
if exceeded a length (10 characters). Along ``XICRA`` there are a few string 
length limitations by different software that need to be sort out from the beginning 
of the process.

See additional details of this name format limitations under user-guide section: :ref:`format fastq files<format-fastq-files>`

Once sample names are checked this module allows the user to copy files into the project 
folder initiate or only link them using a symbolic link to avoid duplicated raw data. 

If using the the ``--project`` option, this module would create a single folder for each
sample and a folder named as ``raw`` including the linked or copied fastq files. See additional
details of project folder organization :ref:`here<project-organization>`.

..
.. TODO: parameters
..

.. include:: ../../links.inc