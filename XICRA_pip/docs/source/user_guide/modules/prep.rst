.. ############################
.. _prep-description:
.. ############################

prep
====

This module prepares fastq files from a sequencing run for later usage.

It initially checks the length of the name and advises the user to rename samples 
if exceeded a length (10 characters). Along ``XICRA`` there are a few string 
length limitations by different software that need to be sort out from the beginning 
of the process.

See additional details of this name format limitations under user-guide section: :ref:`format fastq files<format-fastq-files>`

Once sample names are checked this module allows the user to copy files into the project 
folder initiate or only link them using a symbolic link to avoid duplicated raw data. 


In addition, when multiples files have been generated for the same
sample e.g different lanes this module can merge them. It concatenates these files according the common
identifier and generates a unique file, one per paired-read if necessary, check :ref:`format fastq files<format-fastq-files>`.


If using the the ``--project`` option, this module would create a single folder for each
sample and a folder named as ``raw`` including the linked or copied fastq files. See additional
details of project folder organization :ref:`here<project-organization>`.

.. ##################
.. _run-prep:
.. ##################
How to run the prep module
--------------------------
Executing the following:

.. code-block:: sh

   XICRA prep -h


The different options and parameters for this module should appear in the command line prompt:

.. function:: Module XICRA prep help

   :param -h --help: Show this help message and exit. 
   
.. function:: Module XICRA prep Input/Output

   :param --input: Folder containing the files with reads. Files could be .fastq/.fq/ or fastq.gz/.fq.gz. See -\ -help_format for additional details. REQUIRED.
   :param --output_folder: Output folder. Name for the project folder.
   :param --single_end: Single end files. Default mode is paired-end. Default OFF.
   :param --batch: Provide this option if input is a file containing multiple paths instead a path.      
   :param --in_sample: File containing a list of samples to include (one per line) from input folder(s). Default OFF.
   :param --ex_sample: File containing a list of samples to exclude (one per line) from input folder(s). Default OFF.
   :param --detached: Isolated mode. No project folder initiated for further steps. Default OFF.
   :param --include_lane: Include the lane tag (*L00X*) in the sample identification. See -\ -help_format for additional details. Default OFF.
   :param --include_all: Include all file name characters in the sample identification. See -\ -help_format for additional details. Default OFF.

   :type input: string
   :type output_folder: string
   :type in_sample: string 
   :type ex_sample: string

.. function:: Module XICRA prep options

   :param --threads: Number of CPUs to use. Default: 2.
   :param --copy_reads: Instead of generating symbolic links, copy files into output folder. Default OFF.
   :param --merge_Reads: Merge files corresponding to the same sample. Used in combination with -\ -include_lane and -\ -include_all will produce different results. Please check, -\ -help_format or https://xicra.readthedocs.io/en/latest/user_guide/info/info_index.html
   :param --rename: File containing original name and final name for each sample separated by comma. No need to provide a name for each pair if paired-end files. If provided with option '-\ -merge_Reads', the merged files would be renamed accordingly.
   
   :type threads: int
   :type rename: string
   
.. function:: Module XICRA prep additional information
  
   :param --help_format: Show additional help on name format for files.
   :param --help_project: Show additional help on the project scheme.
   :param --debug: Show additional message for debugging purposes.   
   
- For further information of the module functionallity, check this :doc:`page <../../api/modules/prep>`.


Output of trimm for each sample
-------------------------------
After the ``prep`` module excution, a data folder will be generated. Inside it, for each sample, a new folder will be created 
(called as the sample). These sample folders will contain links to all the raw files that correspond to each sample (if -\ -copy_reads has
been selected, instead of links the copied files). 



.. include:: ../../links.inc