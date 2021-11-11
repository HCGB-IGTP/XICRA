QC
===

This module calls fastqc_, a quality check program, to analyze the 
quality of each sample. 

For each sample, a folder called "fastqc" is generated with the 
quality analysis, reported in an html file.

In addition, by default, a multiQC_ report is generated for all the samples.


Workflow
--------
TODO: build image

Functions
---------

.. automodule:: XICRA.modules.qc
   :members:
   
Other modules
-------------
The ``QC`` module also uses the ``multiQC_report`` module:

.. automodule:: XICRA.scripts.multiQC_report
   :members:

.. include:: ../../links.inc