# XICRA: Small RNAseq pipeline for paired-end reads.

## Table of Contents

- [Decription](#description)
  * [Installation](#installation)
  * [XICRA.stats](#xicra.stats)
  * [Etymology](#etymology)
- [Supplementary information](#supplementary-information)
- [Documentation](#documentation)
- [License](#license)
- [Citation](#citation)

## Description

XICRA is a python pipeline developed in multiple separated modules that it is designed to take paired end fastq reads, trim adapters and low-quality base pairs positions, and merge reads (R1 & R2) that overlap. Using joined reads it describes all major RNA biotypes present in the samples including miRNA and isomiRs, tRNA fragments (tRFs) and piwi associated RNAs (piRNAs). 

So far, XICRA produces a miRNA analysis at the isomiR level using joined reads, multiple software at the user selection and following a standardization procedure. 
Results are generated for each sample analyzed and summarized for all samples in a single expression matrix. This information can be processed at the miRNA or 
isomiR level (single sequence) but also summarizing for each isomiR variant type. This information can be easily accessed using the accompanied R package 
[XICRA.stats](https://github.com/HCGB-IGTP/XICRA.stats). Although the pipeline is designed to take paired-end reads, it also accepts single-end reads. 

The workflow of the pipeline is described here.
![Workflow](workflow/XICRA_pipeline.png "XICRA pipeline")

### Installation

XICRA is available in pip https://pypi.org/project/XICRA/.

To install type:

`pip install XICRA`

### XICRA.stats

We additionally provide a supplementary R package for parsing and plotting some XICRA results. See additional details [here](https://github.com/HCGB-IGTP/XICRA.stats).

Install it in R using:

```R
# Install XICRA.stats version from GitHub:
# install.packages("devtools")
devtools::install_github("HCGB-IGTP/XICRA.stats")
```


### Etymology
XICRA means in Catalan "tassa petita, més aviat alta i estreta, emprada expressament per a prendre la xocolata desfeta o cafè. També és una unitat de mesura de volum per a líquids que es feia servir a Catalunya per a l'oli, vi, o llet. https://ca.wikipedia.org/wiki/Xicra"

## Supplementary Information
In this repository we provide supplementary information for the original paper describing the method. See additional details in folder BMC_bioinformatics_paper or [here](BMC_bioinformatics_paper/README.md)

## Documentation
See a full documentation deatils [here](https://xicra.readthedocs.io/)

## License 
MIT License
Copyright (c) 2020 HCGB-IGTP

See additional details [here](XICRA_pip/LICENSE)

## Citation
Sanchez-Herrero et. al .... 2020


