# Recipe for XICRA simulations

To illustrate the potential of paired-end reads at the isomiR level analysis we have generated computer simulations to test the impact of technical errors from single end or paired-end reads. We followed the guidelines previously described for isomiR computer simulations by [Amsel et al. 2017](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-017-1772-z)

We created biological variation and technical variation using multiple high throughput sequencing profiles and evaluated the performance of the simulation using sensitivity and precision of the isomiRs detected under several circumstances. 

## Biological variation

We created artificial miRNA isoforms from Homo sapiens mature and hairpin sequences in miRBase using the bioinformatic scripts previously described at [isomiR-Benchmark](https://github.com/DanielAmsel/isomiR-Benchmark). We additionally added the canonical fasta sequence of each miRNA to the variant dataset generated. 

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


```

From the variant frequency table generated we select 100 miRNAs and discarded some random variants and generated random distribution frequencies of the variant types generated. To simplify the interpretation and evaluation, for each variant type, we selected a single isoform. For each variant type we included the corresponding frequencies generated to a total amount of 100 sequences for each mature miRNA. 

```sh
## get frequency
python get_freq.py all_seqs.fa fasta all_seqs.freqs.csv
```

## Technical simulation

For each biological dataset generated, we used ART (version Mount Rainier 2016–06-05 ) with Illumina HiSeq2500 and MiSeq-v1 sequencing system profiles using paired-end mode to simulate next generation sequencing (NGS) reads. We grouped all isoforms according to length and generated NGS simulation for each length subset to finally merged them all in a single file for each read accordingly. We used a 10x sequencing coverage for each input fasta sequence. 

As previously noted (Amsel, 2017), due to the nature of the ART simulation we had to parse and omit about half the total reads generated as they were reverse complemented. We only discarded reverse complemented R1 reads and its R2 counterpart accordingly. Due to the implementation based on frequencies that we did in the biological variation procedure, we made sure when applying the same coverage for each sequence, and discarding reverse complement, that the biological variation frequencies generated would be maintained in the NGS simulation. The observed range would vary from 5-500 counts for each single variant type simulated. 

We evaluated the performance of using PE reads for miRNA isomiR analysis using the NGS simulation datasets and the pipeline XICRA. For each dataset, we used the miRNA module using paired-end mode and single end mode for the R1 and R2 reads. For the paired-end mode we initially joined reads using two different join percentage difference cutoff (fastq-join parameter) to test the effect of using 100% perfect R1 and R2 reads or allowing the default difference (8%) along the minimum default overlap length cutoff (6 bp). For the single-end mode, we used the total R1 reads simulated or the total R2, reversed complemented using seqtk software, respectively. We generated a miRNA analysis at the isomiR level using the three different software available within XICRA: miraligner, sRNAbench and optimir.

All these steps mention here are implemented in simulation_sender.py

```sh
## simulate NGS reads using art
python simulation_sender.py --fasta all_seqs.fa --folder example --reads PE --seqSys HS25 --fcov 10 -t 2 --freqs all_seqs.freqs.csv -n 10 -r 10 --art_bin ./art_illumina -m 50 -s 5 --database db_mirbase
```

## Performance evaluation

Using the biological variation frequencies generated as true positives for each dataset, we evaluated the amount of isomiRs detected for each software and each type of read. For paired-end reads, we also used a different percentage difference cutoff. 

For each isomiR, the detected counts were classified as: True positives (TP) when observed counts matched the expected counts; false positives (FP) when observed counts exceeded the expected counts and were wrongly assigned; false negatives (FN) when observed counts did not get to the minimum expected counts. We calculated the sensitivity or recall as TP/(TP+FN) and the precision or specificity as TP/(TP+FP). We also reported True Negatives (TN) when expected counts were not observed and new generation isomiRs when new variants or miRNA appeared and were not expected. 
