.. ############################
.. _biotype-description:
.. ############################

biotype
=======
This module generates a RNA biotype analysis. The aim of this computation 
is to check if there are samples with a very different configuration, outliers
that show different proportions of uniquely mapped, multimapped or no mapped
reads or with different quantities of miRNA, misc_RNA,... than the rest.
If this happens, it could be due to possible differences in sample manipulation, 
extraction, library preparation, sequencing, etc. Those samples should be 
excluded. 

The module is divided in two steps: 

1. Mapping the reads: performed with STAR_ software.

2. Feature counts: perforfed with featureCounts_.

The output of this module is ‘descriptive’, thus, we won’t need the output to continue 
with the analysis, it is an informative step to know the RNA types that are present
in each sample.  

**Note that**: the mapping process performed with STAR software requires high RAM values. 
We are working on the implementation of alternatives to be able to execute the ``biotype``
module in computers with less RAM capacity. 

.. ##################
.. _run-biotype:
.. ##################
How to run the biotype module
------------------------------
Executing the following:

.. code-block:: sh

   XICRA biotype -h


The different options and parameters for this module should appear in the command line prompt:

.. function:: Module XICRA biotype help

    :param -h --help: Show this help message and exit. 
   
.. function:: Module XICRA QC Input/Output

    :param --input: Folder containing the files with reads. Files could be .fastq/.fq/ or fastq.gz/.fq.gz. See -\ -help_format for additional details. REQUIRED.
    :param --output_folder: Output folder. Name for the project folder.
    :param --single_end: Single end files. Default OFF. Default mode is paired-end.
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

.. function:: Module XICRA biotype options

    :param --threads: Number of CPUs to use. Default: 2.
    :param --annotation: Reference genome annotation in GTF format.
    :param --limitRAM: limitRAM parameter for STAR mapping. Default 20 Gbytes.
    :param --noTrim: Use non-trimmed reads [or not containing '_trim' in the name].
    :param --skip_report: Do not report statistics using MultiQC report module. Default OFF. See details in -\ -help_multiqc

   
    :type threads: int
    :type annotation: string
    :type limitRAM: int

.. function:: Module XICRA biotype parameters

    :param --no_multiMapping: Set NO to counting multimapping in the feature count.By default, multimapping reads are allowed. Default: False
    :param --stranded STRANDED: Select if reads are stranded [1], reverse stranded [2] or non-stranded [0], Default: 0.

.. function:: Module XICRA biotype reference genome

    :param --fasta: Reference genome to map reads.
    :param --genomeDir: STAR genomeDir for reference genome.


.. function:: Module XICRA biotype additional information
  
    :param --help_format: Show additional help on name format for files.
    :param --help_project: Show additional help on the project scheme.
    :param --help_RNAbiotype: Show additional help on the RNAbiotype paired-end reads process.
    :param --debug: Show additional message for debugging purposes.   
   
- For further information of the module functionallity, check this :doc:`page <../../api/modules/biotype>`.

Output of biotype
-----------------
Inside the data folder of each sample, a 'map' directory will be generated containing a
report from the mapping of MultiQC. After that, a final report will also be created in the
'report' folder with the featureCounts information of all samples. 

.. include:: ../../links.inc