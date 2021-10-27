######################
Additional information
######################

.. only:: html

    :Version: |version|
    :Date: |today|
    

.. contents::

.. # 
.. TODO: complete this page
..

.. _background:

Project Status
==============

The project is available on `XICRA github`_. You can 
`report issues`_ or contribute to the project by `forking the project`_
and `creating pull requests`_. 
See additional details on XICRA :ref:`developer guidelines<developer>` and how to :ref:`work with Git<git-guidelines>`.


Background
==========

MicroRNAs (miRNAs), a class of small non-coding RNAs (ncRNAs), have an average length of 21–23 nucleotides (nt). 
They have been widely studied as endogenous regulatory molecules that modulate gene expression post-transcriptionally 
by inducing target mRNA silencing and decay. Additional roles beyond negative modulation of mRNA function have 
also been proposed. 

Most miRNA expression studies based on next generation sequencing (NGS) performed to this date have summarized all 
the reads mapping to a specific miRNA locus or miRNA sequence with or without mismatches and assign it to a single 
miRNA entity (a miRBase reference database entry). However, this type of analysis neglects the fact that not all reads 
are identical to the mature reference sequence in miRBase. Small RNA sequencing NGS methodology has revealed that 
miRNAs can frequently appear in the form of multiple sequence variants or isoforms (termed isomiRs).

``XICRA`` is a pipeline to analyze small RNA sequencing (small RNA-seq) data, which allows detecting and quantifying 
isomiR level variation in miRNA.


.. _license:

License
=======

This is brief description of the license terms for ``XICRA``. 


.. _credits:

Credits
=======

Give credit to who deserves


.. _citation:

Citation
========

Please cite the ``XICRA`` project, code or documentation using the following links:

.. github_stats:

Github statistics
=================

These are Github statistics generated for each release of the code.

.. include:: ./github_stats.txt
   :literal:

.. _whats_new:

What's new?
===========

See below information on differences among each release of the code generated.

.. include:: ./whats_new.txt
   :literal:
  
 
.. include:: ../links.inc
    