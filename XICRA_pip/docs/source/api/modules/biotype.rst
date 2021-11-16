biotype
=======
.. 
.. TODO: XICRA.scripts.RNAbiotype functions description
..


This module works in two steps. First, it maps the unjoined reads
to a given reference genome. Secondly, it generates a summary of the RNA
types found in each sample. 

Mapping
-------
The mapping step is performed with STAR_ software and the code is in a separate
script, ``mapReads.py``. To introduce the reference genome the are two options:

1. Introduce the path to the fasta sequence of the human genome, option -\ -fasta, and
   indicate the path to a reference genome annotation in GTF format,  -\ -annotation.
2. Introduce the STAR indexed genome file, option -\ -genomeDir. This directory, 
   called 'STAR_index' is created after the execution of STAR, thus, if it is 
   already generated it can be reused, check the STAR_ documentation.

The -\ -limitRAM parameter is also important. It indicates the maximum RAM (in bytes) that 
the computation will use to prevent the computer collapse. Note that, the STAR software requires high 
values of RAM  in order to do the mapping. Thus, a 'regular' laptop may not be able to perform 
this computation. We are currently working in the implementation of alternatives to this software 
that can be used with any computer. 

As a result, the mapping step generates a 'map' folder with a BAM file and a report from MultiQC
for each sample.

Feature counts
--------------
The second step performs the RNA biotype analysis  using the featureCounts_ software, 
code located in the script ``RNAbiotype.py``. After the computation for each sample, a final
MultiQC report is generated to compare the composition of each sample and, eventually, detect outliyers.


Workflow
--------
TODO: build image 


Functions
---------

.. automodule:: XICRA.modules.biotype
   :members:
   
Other functions
---------------
This module calls other scripts:

.. automodule:: XICRA.scripts.mapReads
   :members:

.. automodule:: XICRA.scripts.RNAbiotype
   :members:

.. include:: ../../links.inc