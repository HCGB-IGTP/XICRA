## download raw data files in raw_data folder
mkdir raw_data

# donwload GSE114923 files
awk '{print $12}' info_project.tsv | tr ';' '\n' | grep -v 'submitted' > ftp_info.txt
wget -i ftp_info.txt

## paired-end XICRA analysis
XICRA prep --input ./raw_data/ --output XICRA_analysis
XICRA trimm --input XICRA_analysis --adapters_a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC --adapters_A GATCGTCGGACTGTAGAACTCTGAACGTGTAGATCTCGGTGGTCGCCGTATCATT --threads 24

## perc_diff_0
XICRA join --input XICRA_analysis/ --perc_diff 0 --threads 24
XICRA miRNA --threads 24 --software miraligner --database db_miRBase --input XICRA_analysis/

## perc_diff_8
XICRA prep --input ./raw_data/ --output XICRA_analysis_perc_diff_8
XICRA trimm --input XICRA_analysis_perc_diff_8 --adapters_a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC --adapters_A GATCGTCGGACTGTAGAACTCTGAACGTGTAGATCTCGGTGGTCGCCGTATCATT --thr$
XICRA join --input XICRA_analysis_perc_diff_8/ --perc_diff 8 --threads 24
XICRA miRNA --threads 24 --software miraligner --database db_miRBase --input XICRA_analysis_perc_diff_8/

## R1
mkdir reads_R1
cd reads_R1
for i in `dir ../XICRA_analysis/data/`; do `ln -s ../XICRA_analysis/data/$i/trimm/*R1*`; done
cd ../
XICRA miRNA --detached --single_end --input reads_R1/ --threads 24 --software miraligner --database db_miRBase/ --output_folder analysis_R1

## R2
mkdir reads_R2
cd reads_R2
for i in `dir ../XICRA_analysis/data/`; do echo "## Reverse complement $i"; R2_name=$i".fastq"; `seqtk seq -r ../XICRA_analysis/data/$i/trimm/*R2* > $R2_name`; done
cd ../
XICRA miRNA --detached --single_end --input reads_R2/ --threads 24 --software miraligner --database db_miRBase/ --output_folder analysis_R2

