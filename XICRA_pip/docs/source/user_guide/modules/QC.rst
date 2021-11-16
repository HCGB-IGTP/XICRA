.. ############################
.. _qc-description:
.. ############################

QC
===

This module calls fastqc_, a quality check program, to analyze the 
quality of each sample. 

By default, creates a final MultiQC_ report with all the samples. It is useful 
to check if there are outliers among the input samples. It can be disabled using
the option -\ -skip_report. 

.. ##################
.. _run-qc:
.. ##################
How to run the QC module
--------------------------
Executing the following:

.. code-block:: sh

   XICRA QC -h


The different options and parameters for this module should appear in the command line prompt:

.. function:: Module XICRA QC help

    :param -h --help: Show this help message and exit. 
   
.. function:: Module XICRA QC Input/Output

    :param --input: Folder containing the files with reads. Files could be .fastq/.fq/ or fastq.gz/.fq.gz. See -\ -help_format for additional details. REQUIRED.
    :param --output_folder: Output folder. Name for the project folder.
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

.. function:: Module XICRA QC options

    :param --skip_report: Do not report statistics using MultiQC report module. Default OFF.
    :param --threads: Number of CPUs to use. Default: 2.
   
    :type threads: int
    :type rename: string

.. function:: Module XICRA QC additional information
  
    :param --help_format: Show additional help on name format for files.
    :param --help_project: Show additional help on the project scheme.
    :param --help_multiqc: Show additional help on the multiQC module.
    :param --debug: Show additional message for debugging purposes.   
   
- For further information of the module functionallity, check this :doc:`page <../../api/modules/QC>`.


Output of QC
------------
After the ``QC`` module excution,in the data folder for each sample, a new directory will be created.
It will be called 'fastq' and will contain the the quality analysis.

In addition, if -\ -skip_report was not selected, a folder called 'report' will also be generated. In
this folder, the different reports that may be created after the execution of other modules (as ``trimm``)
will be stored. These reports will always show information of all the samples. In this case, a subfolder 
called 'fastqc' will be built with the MultiQC report in order to visualize the quality 
analysis of all samples.


.. include:: ../../links.inc