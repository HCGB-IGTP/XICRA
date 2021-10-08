.. ############################
.. _miRNA-description:

miRNA
=====

The ``XICRA`` pipeline provides this module, ``miRNA``, to generate the miRNA analysis.
MicroRNAs (miRNAs), a class of small non-coding RNAs, have an average length 
of 21–23 nucleotides (nt) and modulate gene expression post-transcriptionally. Most miRNA 
expression studies based on next generation sequencing (NGS) data, summarize all the reads 
mapping to a specific miRNA locus or miRNA sequence with or without mismatches and assign
it to a single miRNA entity, this is, a unique miRBase_ reference database entry (miRBase
database is a searchable database of published miRNA sequences and annotation).

However, this type of analysis neglects the fact that **not all reads are identical to the** 
**reference sequence in miRBase, which is called the canonical sequence**.
Small RNA sequencing NGS methodology has revealed that miRNAs can frequently appear in the 
form of multiple sequence variants or isoforms (termed isomiRs). Each isomiR can modulate 
gene expression post-transcriptionally. 


With the ``miRNA`` module we capture all the available information at different levels, 
the analysis can be done at miRNA or isomiR level. 

Its functionalyty is divided in three steps: 
   #. miRNA analysis
   #. Standarization of the results 
   #. Expression count matrix generation with miRTop_ (Command line tool to annotate with a standard naming miRNAs e isomiRs).


The analysis can be performed with three different softwares:

- Miraligner_: maps small RNA data to miRBase repository.
- Optimir_: algorithm for integrating available genome-wide genotype data into miRNA sequence alignment analysis.
- sRNAbench_: application for processing small-RNA data obtained from NGS platforms. 
  Unfortunately, the downloading of this tool is no longer available, thus, only users 
  with ``sRNAbench`` already installed will be able to run the ``XICRA`` analysis with it. 

According to our tests, published in this article_, ``miraligner`` is the software with 
the best performance for the miRNA analysis. 

How to run the miRNA module
---------------------------
Executing the following:

.. code-block:: sh

   XICRA miRNA -h


The different options and parameters for this module should appear in the command line prompt:

.. function:: Module XICRA miRNA optional arguments
      
   :param -h --help: Show this help message and exit. 


.. function:: Module XICRA miRNA Input/Output

   :param --input: Folder containing a project or reads, according to the mode selected. Files could be .fastq/.fq/ or fastq.gz/.fq.gz. See --help_format for additional details. REQUIRED.
   :param --output_folder: Output folder.
   :param --single_end: Single end files. Default mode is paired-end. Default OFF.
   :param --batch: Provide this option if input is a file containing multiple paths instead a path.      
   :param --in_sample: File containing a list of samples to include (one per line) from input folder(s). Default OFF.
   :param --ex_sample: File containing a list of samples to exclude (one per line) from input folder(s). Default OFF.
   :param --detached: Isolated mode. --input is a folder containing fastq reads. Provide a unique path o several using --batch option.
   :param --include_lane: Include the lane tag (*L00X*) in the sample name. See --help_format for additional details. Default OFF.
   :param --include_all: Include all characters as tag name before read pair, if any. See --help_format for additional details. Default OFF.
   :param --noTrimm: Use non-trimmed reads (or not containing '_trim' in the name).
   
   :type input: string
   :type output_folder: string
   :type in_sample: string 
   :type ex_sample: string

      
.. function:: Module XICRA miRNA options

   :param --threads: Number of CPUs to use. Default: 2. 
   :param --species: Species tag ID. Default: hsa (Homo sapiens).
   :param --database: Path to store miRNA annotation files downloaded: miRBase, miRCarta, etc.
   :param --miRNA_gff: miRBase hsa GFF file containing miRNA information.
   :param --hairpinFasta: miRNA hairpin fasta file.
   :param --matureFasta: miRNA mature fasta file.
   :param --miRBase_str: miRBase str information.
   
   :type threads: int 
   :type species: string 
   :type database: string
   :type miRNA_gff: string
   :type hairpinFasta: string
   :type matureFasta: string
   :type miRBase_str: string
   
   
.. function:: Module XICRA miRNA software

   :param --software: Software to analyze miRNAs, sRNAbench, optimir, miraligner. Provide several input if desired separated by a space. REQUIRED.
   
   :type software: string

.. function:: Module XICRA miRNA additional information
  
   :param --help_format: Show additional help on name format for files.
   :param --help_project: Show additional help on the project scheme.
   :param --help_miRNA: Show additional help on the miRNA paired-end reads process.
   :param --debug: Show additional message for debugging purposes.


Output of miRNA for each sample
-------------------------------
As the rest of the modules, the ``miRNA`` module will generate a folder in each of the 
sample directories called "miRNA". Inside this folder another two will be created for 
each of the softwares selected. For example, if we have executed ``--software optimir miraligner``
we will obtain four output folders:

- data/sampleName/miRNA/optimiR
- data/sampleName/miRNA/miraligner 
- data/sampleName/miRNA/optimiR_miRTop
- data/sampleName/miRNA/miraligner_miRTop

The first two folders will store the outputs of the corresponding softwares in their particular
format.

The folders ended in "_miRTop" will contain the results in the miRTop standarized format. 

Finally, the expression count matrix will be stored in .tsv format. Following the previous 
example, these files would be located in:

- data/sampleName/miRNA/optimiR_miRTop/counts/mirtop.tsv
- data/sampleName/miRNA/miraligner_miRTop/counts/mirtop.tsv


Expression count matrix for each sample
---------------------------------------

As a result for each sample (and software used) we will end up with a table like this, mirtop.tsv, 
called the expression count matrix. Here we can see an example of this matrix:

..
.. TODO: resize the table
..

.. csv-table::
   :file: mirtop_example.csv
   :header-rows: 1

- UID: unique identifier (UID) for each sequence defined by miRTop.
- Read: DNA sequence.
- miRNA: miRNA precursor, identifier defined by miRBase for each miRNA canonical sequence.
- Variant: variant type of each isomiR, 'NA' for the canonical sequence (checkout the miRTop `variant nomenclature`_).
- The following four columns indicate the amount of base pairs added or substracted, compared to the canonical sequence.
- SampleName: raw read count expression for this sample.

Output of miRNA, comparing samples
----------------------------------

On the other hand, as other modules, ``miRNA`` also builds an output to compare samples. In the folder
report/miRNA, three different files will be created for each software executed. For example, if we have run 
``--software miraligner``, we will obtain the following files: 

- **report/miRNA/miRNA_expression-miraligner_dup.csv**: Matrix with the number of reads of each UID of each sample that are duplicated. 
  Normally, they occur when some bases are added at the beginning and the end, so it cannot be differentiated if 
  they are 3p:+1;5p:+2 or 3p:+2;5p:+1. In those cases, they will both have the same UID. They are removed 
  (they typically have very few counts).
- **report/miRNA/miRNA_expression-miraligner.csv**: Final matrix (without the duplicated UIDs). Number of counts of each UID of each sample,
  to be further analyzed with ``R``. 
- **report/miRNA/miRNA_expression-miraligner_seq.csv**: table with the DNA sequence corresponding to each UID. 

The analysis of the matrix stored in miRNA_expression-miraligner.csv can be done at the isomiR level, differenciating by
UID, variant type or miRNA (just considering the miRNA identifier).  It can be done with the package XICRA.stats_.

.. include:: ../../links.inc