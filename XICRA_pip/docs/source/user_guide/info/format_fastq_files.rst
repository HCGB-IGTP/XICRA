.. ########################
.. _format-fastq-files:
.. ########################

Sample name details
===================

``XICRA`` accepts a wide range of files names. 

Format for fastq files can be:

- name.fastq.gz
- name_1.fastq.gz, adding '1' or '2' to specify the read
- name_R2.fastq.gz, adding 'R1' or 'R2' to specify the read
- name_L001_R1.fastq.gz, adding the lane information as L00x after the name
- name_L001_R1_001.fastq.gz, adding 00X at the end. This naming is useful when the fastqfiles of the sample sample had been cut in different files).
- name_L001_XYZ_R1_001.fastq.gz, there can be extra info for each file.

The file names should be structures as follows:

Length limitation
-----------------

There is a limitation for the sample ID ('name') of 25 characters.
** XICRA provides an option to rename samples if necessary: module prep option --rename. **

Extensions
----------
The suported extensions are:

- name_L00x_R2.fastq   
- name_L00x_R2.fq
- name_L00x_R2.fastq.gz   
- name_L00x_R2.fq.gz

Single end files
----------------

It is possible to provide NGS single-end files although some steps of the process could be accomplished using single-end files.

- name.fastq.gz
- name.fastq
- name.fq

** Use option --single-end in the different XICRA modules. **

Paired-end files
----------------

Paired-end files are full supported. The format for these files are:

- name_1.fastq.gz or name_R1.fastq.gz
- name_2.fastq.gz or name_R2.fastq.gz

** No parameter is needed in this case.** 

Lane information
----------------

In some cases, files might contain lane information (*L00x* and/or *00x*).
``XICRA`` supports these names as long as follow these examples:

- name_L00x_R1.fastq.gz   name_L00x_R2.fastq.gz
- name_L00x_1.fastq.gz    name_L00x_2.fastq.gz


Name extensions
---------------

It can also be the case that the reads of a sample are divided in different files. In those cases, the files
should contain a name final extension: 

- sample1_L001_R1_001.fastq.gz sample1_L001_R2_001.fastq.gz
- sample1_L001_R1_002.fastq.gz sample1_L001_R2_002.fastq.gz
- sample1_L002_R1_001.fastq.gz sample1_L002_R2_001.fastq.gz
- sample1_L002_R1_002.fastq.gz sample1_L002_R2_002.fastq.gz

Extra information
-----------------
In some cases, files might contain other extra information. In the following example, XYZ is the extra information:

- sample1_L001_XYZ_R1_001.fastq.gz sample1_L001_XYZ_R2_001.fastq.gz
- sample1_L001_XYZ_R1_002.fastq.gz sample1_L001_XYZ_R2_002.fastq.gz


Sample identification
=====================

``XICRA`` will store the names of all the input files. After that, it will identify the samples. 
It can be the case that more than one file belong to the same sample. In order to pass this information to 
``XICRA`` a combination of the following parameters may be needed depending on the input file name characteristics:

--include_lane
--------------

If you want to include lane tags (*L00X*) into each  each sample name (differentiate samples considering the lane):
** Use option --include-lane within each module and the lane tag will also be used to identify samples **

However, if you want to consider as a single sample the different lanes, you need to merge the fastq files from the 
different lanes:
** Use option --merge_Reads within module prep **

As an example, considering the input files:

- sample1_L001_R1.fastq.gz   sample1_L001_R2.fastq.gz
- sample1_L002_R1.fastq.gz   sample1_L002_R2.fastq.gz
- sample1_L003_R1.fastq.gz   sample1_L003_R2.fastq.gz
- sample1_L004_R1.fastq.gz   sample1_L004_R2.fastq.gz

   #. By adding the option ``--include_lane`` within module ``prep``:
   
      - sample1_L001_R1.fastq.gz   sample1_L001_R2.fastq.gz
      - sample1_L002_R1.fastq.gz   sample1_L002_R2.fastq.gz
      - sample1_L003_R1.fastq.gz   sample1_L003_R2.fastq.gz
      - sample1_L004_R1.fastq.gz   sample1_L004_R2.fastq.g
      ``XICRA`` will identify four samples.
      ** Remember to use option --include_all within each module**
      
   #. By adding the options ``--include_lane --merge_Reads`` within module ``prep``:
   
      - sample1_R1.fastq.gz  sample1_R2.fastq.gz
      ``XICRA`` will only identify one sample.
      

--include_all
-------------

In some cases, files might contain other extra information and it is necessary to include all the information of the 
file name to identify samples:
** In that case use ``--include-all``.** 

If you need to merge fastq files that only differ in the final extension (_001, _002, ...), 
use option ``--merge_Reads`` within module prep together with ``--include-all`` **


As an example, considering the input files:

- sample1_L001_XYZ_R1_001.fastq.gz sample1_L001_XYZ_R2_001.fastq.gz
- sample1_L001_XYZ_R1_002.fastq.gz sample1_L001_XYZ_R2_002.fastq.gz
- sample1_L002_XYZ_R1_001.fastq.gz sample1_L002_XYZ_R2_001.fastq.gz
- sample1_L002_XYZ_R1_002.fastq.gz sample1_L002_XYZ_R2_002.fastq.gz

   #. By adding the options ``--include_all --merge_Reads`` within module ``prep``:
   
      - sample1_L001_XYZ_R1.fastq.gz  sample1_L001_XYZ_R2.fastq.gz
      - sample1_L002_XYZ_R1.fastq.gz  sample1_L002_XYZ_R2.fastq.gz
      ``XICRA`` will identify two samples
      ** Remember to use option --include_all within each module**
      




