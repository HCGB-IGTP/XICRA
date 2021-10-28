.. ########################
.. _format-fastq-files:
.. ########################

Sample name details
===================

``XICRA`` accepts a wide range of file names. The format for fastq files can be:

- name.fastq.gz
- name_1.fastq.gz, with '1' or '2' to specify the read
- name_R2.fastq.gz, 'R1' or 'R2' to specify the read
- name_L001_R1.fastq.gz, with the lane information as 'L00X' or '00X' after the name
- name_L001_R1_001.fastq.gz, with '00X' at the end of the file name. This naming is used 
  when the fastq of a sample had been cut in different files.
- name_L001_XYZ_R1_001.fastq.gz, there can be extra info for each file, XYZ.

The input file names should be structured considering the following aspects:

Length limitation
-----------------

There is a limitation for the sample ID ('name') of 25 characters.

``XICRA`` provides an option to rename samples if necessary in the module ``prep``, option **-**\ **-rename**.

Extensions
----------
The suported extensions are:

- name_L00x_R2.fastq   
- name_L00x_R2.fq
- name_L00x_R2.fastq.gz   
- name_L00x_R2.fq.gz

Single-end files
----------------

It is possible to provide NGS single-end files although some steps of the process could not be accomplished 
using single-end files. 

- name.fastq.gz
- name.fastq
- name.fq

Use option **-**\ **-single-end in the different XICRA modules.**

Paired-end files
----------------

Paired-end files are full supported. The format for these files are:

- name_1.fastq.gz or name_R1.fastq.gz
- name_2.fastq.gz or name_R2.fastq.gz

No parameter is needed in to specify this kind of files.

Lane information
----------------

Files might contain lane information (*L00x* and/or *00x*).
``XICRA`` supports these names as long as follow these examples:

- name_L00x_R1.fastq.gz, name_L00x_R2.fastq.gz
- name_L00x_1.fastq.gz, name_L00x_2.fastq.gz


Name extensions
---------------

It can also be the case that the reads of a sample are divided in different files. In those cases, the files
should contain a name final extension: 

- name1_L001_R1_001.fastq.gz, name1_L001_R2_001.fastq.gz
- name1_L001_R1_002.fastq.gz, name1_L001_R2_002.fastq.gz
- name1_L002_R1_001.fastq.gz, name1_L002_R2_001.fastq.gz
- name1_L002_R1_002.fastq.gz, name1_L002_R2_002.fastq.gz

Extra information
-----------------
In some cases, files might contain other extra information. In the following example, XYZ is the extra information:

- name1_L001_XYZ_R1_001.fastq.gz, name1_L001_XYZ_R2_001.fastq.gz
- name1_L001_XYZ_R1_002.fastq.gz, name1_L001_XYZ_R2_002.fastq.gz


Sample identification
=====================

``XICRA`` will store the names of all the input files. After that, it will identify the samples. 
It can be the case that more than one file belong to the same sample. In order to pass this information to 
``XICRA``, a combination of the following parameters may be needed depending on the characteristics of the
input file names:

Option: include_lane
--------------------

If you want to include lane tags (*L00X*, *00X*) into each  each sample name (differentiate samples considering the lane):
Use option **-**\ **-include_lane within each module** and the lane tag will also be used to identify samples.

However, if you want to consider as a single sample the different lanes, you need to merge the corresponding fastq files:
Use option **-**\ **-merge_Reads** within module ``prep``.

As an example, considering the input files:

- name1_L001_R1.fastq.gz, name1_L001_R2.fastq.gz
- name1_L002_R1.fastq.gz, name1_L002_R2.fastq.gz
- name1_L003_R1.fastq.gz, name1_L003_R2.fastq.gz
- name1_L004_R1.fastq.gz, name1_L004_R2.fastq.gz

   #. By adding the option ``--include_lane`` in **all modules**, ``XICRA`` will identify four samples:
   
      - Sample 1: name1_L001_R1, name1_L001_R2
      - Sample 2: name1_L002_R1, name1_L002_R2
      - Sample 3: name1_L003_R1, name1_L003_R2
      - Sample 4: name1_L004_R1, name1_L004_R2
      
      Remember to use option **-**\ **-include_lane within each module**.
      
   #. By adding the options ``--include_lane --merge_Reads`` **within module prep**, ``XICRA`` will only 
      identify one sample, merging all the corresponding files:
   
      - Sample 1: sample1_R1, sample1_R2
      
      

Option: include_all
-------------------

In some cases, files might contain other extra information and it is necessary to use all the information of the 
file name to identify samples:

If that is the case use **-**\ **-include_all in al modules** .

If you want to merge fastq files that only differ in the final extension (_001, _002, ...):

Use options **-**\ **-merge_Reads** **-**\ **-include_all within module prep** and only 
**-**\ **-include_all in the rest of the modules**.

As an example, considering the input files:

- name1_L001_XYZ_R1_001.fastq.gz, name1_L001_XYZ_R2_001.fastq.gz
- name1_L001_XYZ_R1_002.fastq.gz, name1_L001_XYZ_R2_002.fastq.gz
- name1_L002_XYZ_R1_001.fastq.gz, name1_L002_XYZ_R2_001.fastq.gz
- name1_L002_XYZ_R1_002.fastq.gz, name1_L002_XYZ_R2_002.fastq.gz

   #. By adding the option ``--include_all`` in **all modules**, ``XICRA`` will identify four samples:
   
      - Sample 1: name1_L001_XYZ_R1_001, name1_L001_XYZ_R2_001
      - Sample 2: name1_L001_XYZ_R1_002, name1_L001_XYZ_R2_002
      - Sample 3: name1_L002_XYZ_R1_001, name1_L002_XYZ_R2_001
      - Sample 4: name1_L002_XYZ_R1_002, name1_L002_XYZ_R2_002
      
      Remember to use option **-**\ **-include_all within each module**.

   #. By adding the options ``--include_all --merge_Reads`` **within module prep**, ``XICRA`` will identify two samples:
   
      - Sample 1: name1_L001_XYZ_R1, name1_L001_XYZ_R2
      - Sample 2: name1_L002_XYZ_R1, name1_L002_XYZ_R2
      
      Remember to use option **-**\ **-include_all within each module**.
      




