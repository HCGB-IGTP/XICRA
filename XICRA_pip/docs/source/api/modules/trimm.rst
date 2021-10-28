trimm
=====

This module cuts the introduced adapter sequences of the input reads.

For each sample, a folder called "trimm" is generated with the trimmed
reads, called as the original ones +"_trimm".

In addition, a multiQC report is generated.


Workflow
--------
TODO: build image


Functions
---------

.. automodule:: XICRA.modules.trimm
   :members:
   
Other modules
-------------
The ``trimm`` module also uses the ``multiQC_report`` module:

.. automodule:: XICRA.scripts.multiQC_report
   :members:

.. include:: ../../links.inc