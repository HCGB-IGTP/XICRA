# Recipe for XICRA simulations

## Simulated isomiR categories (isomiR-Benchmark)
Five add: FA
Five del: FS
Non-template: NT
SNP Rest: SR
SNP Seed: SS
Three add: TA
Three del: TS

### Variant types for miRTop: 

iso_5p/iso_3p:+/-N. 
 (+) indicates the start is shifted to the right. 
 (-) indicates the start is shifted to the left. 
N the number of nucleotides of difference. 
For instance, if the sequence starts 2 nts after the reference miRNA, the label will be: iso_5p:+2, but if it starts before, 
the label will be iso_5p:-2.

iso_add3p:N. Number of non-template nucleotides added at 3p.
iso_add5p:N. Number of non-template nucleotides added at 5p.

iso_snv_seed: when affected nucleotides are between [2-7].
iso_snv_central_offset: when affected nucleotides is at position [8].
iso_snv_central: when affected nucleotides are between [9-12].
iso_snv_central_supp: when affected nucleotides are between [13-17].
iso_snv: anything else.


## Conversion:
Five add: FA iso_5p:-1
Three add: TA iso_3p:+1

Five del: FS iso_5p:+1
Three del: TS iso_3p:-1

Non-template: NT (iso_add5p / iso_add3p)
SNP Rest: SR: iso_snv_central & iso_snv_central_supp
SNP Seed: SS: iso_snv_seed & iso_snv_central_offset




```sh
mkdir mirBase
cd mirBase

## download latest
wget -nd ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.gz
wget -nd ftp://mirbase.org/pub/mirbase/CURRENT/hairpin.fa.gz

# extract gunzip files
gunzip *

## get hsapiens mature and hairpin
grep 'sapiens' mature.fa | awk '{print $1}' | sed 's/>//' > mature_hsapiens.ids.txt
grep 'sapiens' hairpin.fa | awk '{print $1}' | sed 's/>//' > hairpin_hsapiens.ids.txt

## get hsapiens fasta files
seqtk subseq hairpin.fa hairpin_hsapiens.ids.txt > hairpin_hsapiens.fa
seqtk subseq mature.fa mature_hsapiens.ids.txt > mature_hsapiens.fa

## Simulate isomiRs and reads
cd ../
mkdir isomiR_simulations
cd isomiR_simulations

## simulate isomiR
perl ~/isomiR-Benchmark/create_isomiRs.pl ../mirBase/mature_hsapiens.fa ../mirBase/hairpin_hsapiens.fa

## Concat all sequences & change character
cd ..
cat isomiR_simulation/*fa | sed 's/|/::/' > all_seqs.fa

# add canonical & conver U->T & rename/add character
python rename_canonical.py mature_hsapiens.fa mature_hsapiens_renamed.fa
cat mature_hsapiens_renamed.fa >> all_seqs.fa

## get frequency
python get_freq.py all_seqs.fa fasta all_seqs.freqs.csv

## simulate NGS reads using art
python simulation_sender.py --fasta all_seqs.fa --folder example --reads PE --seqSys HS25 --fcov 10 -t 2 --freqs all_seqs.freqs.csv -n 10 -r 10 --art_bin ./art_illumina -m 50 -s 5 --database db_mirbase
```	

