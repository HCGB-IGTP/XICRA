# XICRA: Small RNAseq pipeline for paired-end reads

## Description

XICRA is a python pipeline developed in multiple separated modules that it is designed to take 
paired end fastq reads, trim adapters and low-quality base pairs positions, and merge reads (R1 & R2) 
that overlap. Using joined reads it describes all major RNA biotypes present in the samples including 
miRNA and isomiRs, tRNA fragments (tRFs) and piwi associated RNAs (piRNAs).

So far, XICRA produces a miRNA analysis at the isomiR level using joined reads, multiple software at the 
user selection and following a standardization procedure. Results are generated for each sample analyzed and 
summarized for all samples in a single expression matrix. This information can be processed at the miRNA or 
isomiR level (single sequence) but also summarizing for each isomiR variant type. This information can be 
easily accessed using the accompanied R package XICRA.stats. Although the pipeline is designed to take 
paired-end reads, it also accepts single-end reads.

## Installation

XICRA is available in the Python Package Index ([pip](https://pypi.org/project/XICRA/)).

To install type: 

`pip install XICRA`

## Documentation

See a full documentation, user guide and manual in [here](https://readthedocs.org/)

## Example
Here we include a brief example on how to use XICRA.

First, we create a python environment and will install XICRA and dependencies.


```sh
## create enviroment
python3 -m venv XICRA_env

## activate it
source XICRA_env/bin/activate

## install XICRA and dependencies
pip install XICRA

## execute XICRA
XICRA -h
```

Now, we can test XICRA by using an example of 100 miRNA simulated and provideded within the repository as an example of simulation.

```sh
## run XICRA example
ln -s ~/BMC_bioinformatics_paper/simulation/example/reads/

## prepare reads
XICRA prep --input reads/ --output_folder test_XICRA

## join reads
XICRA join --input test_XICRA --noTrim

## create miRNA analysis
XICRA miRNA --input test_XICRA --software miraligner sRNAbench

## explore results
ls test_XICRA/report/
```

## License 
MIT License
Copyright (c) 2020 HCGB-IGTP

See additional details [here](LICENSE)

## Citation
Sanchez-Herrero et. al .... 2020

## Authors
Antonio Luna de Haro (v0.1)

Jose F Sanchez-Herrero (v1.0)	
