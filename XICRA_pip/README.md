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

XICRA will require python v3.7 and java (we tested in openjdk 14 2020-03-17).

The XICRA python pipeline is available in `pip` and also available using `conda`.

XICRA depends on multiple third party software that we have listed below.

### Dependencies 

Python XICRA module will install itself along some python modules dependencies (pandas, multiqc, pybedtools, biopython etc.). 

But additionally, XICRA depends on third party software that we listed in the following [table](https://github.com/HCGB-IGTP/XICRA/blob/master/XICRA_pip/XICRA/config/software/soft_dependencies.csv).

### Conda environment

We encourage you to install XICRA and all dependencies using the `conda` environment we created and following these instructions. 

To create a new conda environment, install third party software, install XICRA and missing dependencies, do as follows: 

1) Get requirements file from XICRA git repo

```sh
wget https://raw.githubusercontent.com/HCGB-IGTP/XICRA/master/XICRA_pip/devel/conda/environment.yml
```

2) Create environment and install required packages using conda: 

```sh
conda env create -n XICRA -f environment.yml
```

3) Activate environment and install XICRA
```sh
## activate
conda activate XICRA

## install latest python code
pip install XICRA
```

4) Install missing software:  Unfortunately, a couple of executables are not available neither as a `conda` or `pip` packages. These packages are `miraligner` and `sRNAbench`. We have generated a `shell` script to retrieve and include within your `conda environment`.

```sh
## install missing software
sh XICRA_pip/XICRA/config/software/installer.sh
```

To check everything is fine, try executing the `config` module:
```sh
XICRA config
```

### Python environment

If you are not using a `conda` environment as you might have previously installed all dependencies, we encourage you to create a python environment containing all python modules required for XICRA. See as an example this code:

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

## Documentation

See a full documentation, user guide and manual in [here](https://readthedocs.org/)

## Example
Here we include a brief example on how to use XICRA.

First, we create a python environment and will install XICRA and dependencies. See example details shown before.
Then, we can test XICRA by using an example of 100 miRNA simulated and provideded within the repository as an example of simulation.

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

## Documentation
For a full documentation and details visit Read the Docs site [here](https://xicra.readthedocs.io/). 

See a brief example on how to install and run XICRA [here](https://github.com/HCGB-IGTP/XICRA/tree/master/XICRA_pip#example)

## License 

MIT License

Copyright (c) 2020 HCGB-IGTP

See additional details [here](XICRA_pip/LICENSE)

Developed and maintained by Jose F. Sanchez-Herrero and Lauro Sumoy at HCGB-IGTP

http://www.germanstrias.org/technology-services/genomica-bioinformatica/

## Citation
Sanchez Herrero, J.F., Pluvinet, R., Luna de Haro, A. et al. Paired-end small RNA sequencing reveals a possible overestimation in the isomiR sequence repertoire previously reported from conventional single read data analysis. BMC Bioinformatics 22, 215 (2021). https://doi.org/10.1186/s12859-021-04128-1

## Authors
Antonio Luna de Haro (v0.1)
Jose F Sanchez-Herrero (v1.0)	
