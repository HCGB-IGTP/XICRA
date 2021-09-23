miRNA
=====

..
.. TODO: link trimm module, link miRNA from user guide
..

This module analyses the joined or single reads, that have been 
previously trimmed of different samples. To do so, it uses three 
possible softwares, which can be run at the same time if desired. 

For each sample, an expression matrix using the miRTop nomenclature
is generated, containing the information of the counts at miRNA or
isomiR level; also describing the variant type of each isomiR.

In addition, an expression matrix comparing all the samples is created.


Workflow
--------
TODO: build image


Functions
---------

.. automodule:: XICRA.modules.miRNA
   :members:
   
Other functions
---------
This module calls the function:

.. currentmodule:: XICRA.scripts.generate_DE
.. autofunction:: generate_DE

.. include:: ../../links.inc