.. ########################
.. _format-fastq-files:
.. ########################

Sample name details
===================

``XICRA`` accepts a wide range of files names. 

Format for fastq files can be: name.fastq.gz, name_1.fastq.gz, name_R2.fastq.gz, name_L001_R1.fastq.gz, etc.

Here we provide some guidelines of the format.


Length limitation
-----------------

There is a limitation for the sample ID ('name') of 10 characters.
** XICRA provides an option to rename samples if necessary: module prep option --rename **

Single end files
----------------

It is possible to provide NGS single-end files although some steps of the process could be accomplished using single-end files.
** Use option --single-end in the different XICRA modules **

name.fastq.gz

name.fastq

name.fq

Paired-end files
----------------

Paired-end files are full supported. The format for these files are:

name_1.fastq.gz or name_R1.fastq.gz

name_2.fastq.gz or name_R2.fastq.gz


Lane information:
-----------------

In some cases, files might contain lane information (*L00x* and/or *00x*).
``XICRA`` supports these names as long as follow these examples:

name_L00x_R1.fastq.gz   name_L00x_R2.fastq.gz

name_L00x_1.fastq.gz name_L00x_2.fastq.gz

name_L00x_R1_00x.fastq.gz  name_L00x_R2_00x.fastq.gz

** If you need to merge fastq files from different lanes, use option --merge within module prep **


Extensions:
-----------

name_L00x_R2.fastq   

name_L00x_R2.fq

name_L00x_R2.fastq.gz   

name_L00x_R2.fq.gz


