.. ################################################
.. _quickStart:

Quick start
***********

This is quick start guide for impatient people. 

Installation
============
Build and activate a ``conda`` environment and get ``XICRA`` repository:

.. code-block:: sh

   ## clone repo
   git clone https://github.com/HCGB-IGTP/XICRA.git
   
   ## create conda environemt
   conda env create -n XICRA -f XICRA_pip/devel/conda/environment.yml
   
   ## activate the environment 
   conda activate XICRA
   
   ## install latest python code
   pip install XICRA
   
   ## install missing software
   sh XICRA_pip/XICRA/config/software/installer.sh
   
Sometimes, the activation of the environment is done with ``source activate XICRA`` 
instead of ``conda activate XICRA``. 
Check everything is fine by executing the ``config`` module:

.. code-block:: sh

   XICRA config
   

Execute XICRA
=============
The functionallity of ``XICRA`` is divided in modules. Each of them in charge of 
different parts of the analysis. Here we show how to run each of the
modules:

   #. **Prepation of the data**: ``prep``, ``QC``, ``trimm``, ``join``.
   
      - ``XICRA prep --input input_folder --output output_folder``: preparation of the data.
      - ``XICRA QC --input input_folder``: quality analysis of the samples.
      - ``XICRA trimm --input input_folder --adapters_a XXXXX --adapters_A YYYYYYYY``: trimming the adapters. 
         * ``--adapters_a XXXXX``: sequence of the 3' adapter of R1
         * ``--adapters_A  YYYYYYYY``: sequence of the 3' adapters of R2.
      - ``XICRA join --input input_folder``: join paired end reads.

   #. **Quantification of RNA types**: ``biotype``
   
      - ``XICRA biotype --input input_folder --threads X --fasta file.fa --annotation file.gtf --limitRAM YYYY`` 
         * ``--threads X``: number of cores enabled to run the analysis.
         * ``--fasta file.fa``: fasta file with the reference genome to map reads.
         * ``--annotation file.gtf``: reference genome annotation in GTF format. 
         * ``--limitRAM YYYY``: given RAM in bytes to run the analysis. **Note that:** this process consumes high RAM values.

   #. **miRNA analysis**: ``miRNA``
   
      - ``XICRA miRNA --input input_folder --threads X --software analyisis_tool``
         * ``--threads X``: number of cores enabled to run the analysis.
         * ``--software analyisis_tool``: three options available, ``miraligner``, ``optimir`` and ``sRNAbench``.    

Take into account that these are basic examples, each of the modules can be run with other 
different parameters. To see all the possible parameters of each module run ``-h``. For 
example, for the ``trimm`` module:

.. code-block:: sh

   XICRA trimm -h
   

Test example
------------

Here we include a brief example on how to use ``XICRA``.

First, we create a ``conda`` environment and install ``XICRA`` and its dependencies. 
See example details shown before. Then, we can test ``XICRA`` by using an example 
of 100 miRNA simulated and provided within the repository as an example of simulation.

.. code-block:: sh

   ## run XICRA example
   ln -s ~/BMC_bioinformatics_paper/simulation/example/reads/

   ## prepare reads
   XICRA prep --input reads/ --output_folder test_XICRA

   ## join reads
   XICRA join --input test_XICRA --noTrim

   ## create miRNA analysis
   XICRA miRNA --input test_XICRA --software miraligner optimir

   ## explore results
   ls test_XICRA/report/
   
As a result, we will end up with a folder for each of the run analysis for every sample. 
Thus, in the ``miRNA`` folder, we will obtain the miRNA results for our
samples, performed with two different softwares ``miraligner`` and ``optimir``. 


User data
---------

In the presented example, nor ``QC``, neither ``trimm`` steps were necessary (that is why 
``--noTrim`` was added in the ``join`` command). However, with real data, running ``QC`` 
is highly recommended to check the quality of the samples (and filter outliers if necessary).
Running the ``trimm`` command to eliminate the adapters will be necesary for NGS data. 

On the other hand, the ``biotype`` was skipped as well. This step is only informative: it maps 
and annotates the reads and quantifies the RNA types present in each of the samples. 
**Note that:** this step is very time and RAM consuming.

The ``miRNA`` analysis can be performed whith any of the three available softwares ``miraligner``,
``optimir`` and ``sRNAbench`` (they can be combined as seen in the example). Unfortunately, the
donwloading of ``sRNAbench`` is no longer available, thus, only users with the software
already installed will be able to run the miRNA analysis with it. 

Finally, ``XICRA`` is also able to analyse single end reads, in this case ``--single_end`` should
be added in the commands (and no ``join`` step would be necessary). 

.. ###########
.. _deactivate-env:
.. ###########

Deactivate environment
----------------------

After finished the execution of any ``XICRA`` module or script, it is convenient 
to deactivate the environment. You can just close the terminal but it would be more appopiate
to conveniently deactivate the environment first.

To do so, execute one of the following commands:

.. code-block:: sh
  
   conda deactivate XICRA 
   
.. code-block:: sh
  
   source deactivate XICRA 
