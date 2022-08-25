echo "`date` ... Starting ..."
echo ""

mkdir logs

##########################################
## test single end
##########################################
echo "# ------------------------------ #"
echo "XICRA prep -i ./subset_SE/ -o XICRA_analysis --single_end"
echo "..."
echo ""
XICRA prep -i ./subset_SE/ -o XICRA_analysis --single_end | tee logs/XICRA_analysis.prep.log
echo ""

echo "# ------------------------------ #"
echo "XICRA QC -i XICRA_analysis --single_end --threads 4"
echo "..."
XICRA QC -i XICRA_analysis --single_end --threads 4 | tee logs/XICRA_analysis.qc.log
echo ""

echo "# ------------------------------ #"
echo "XICRA trim -i XICRA_analysis --single_end --threads 4 --adapters_a TGGAATTCTCGGGTGCCAAGG"
echo "..."
XICRA trim -i XICRA_analysis --single_end --threads 4 --adapters_a TGGAATTCTCGGGTGCCAAGG | tee logs/XICRA_analysis.trim.log
echo ""

echo "# ------------------------------ #"
echo "XICRA miRNA -i XICRA_analysis --single_end --threads 4 --software miraligner"
echo "..."
XICRA miRNA -i XICRA_analysis --single_end --threads 4 --software miraligner | tee logs/XICRA_analysis.miRNA.log
echo ""
echo "# ------------------------------ #"
echo ""
echo ""
echo ""   
echo "`date` ... Single End Test finished ..."
echo ""
echo ""
echo ""
echo ""
##########################################

##########################################
## test paired-end
##########################################

echo "`date` ... Starting PE test..."

echo "# ------------------------------ #"
echo "XICRA prep -i ./subset_PE/ -o XICRA_analysis_PE"
echo "..."
XICRA prep -i ./subset_PE/ -o XICRA_analysis_PE | tee logs/XICRA_analysis_PE.prep.log
echo ""

echo "# ------------------------------ #"
echo "XICRA QC -i XICRA_analysis_PE --threads 4"
echo "..."
XICRA QC -i XICRA_analysis_PE --threads 4 | tee logs/XICRA_analysis_PE.qc.log
echo ""

# ATTENTION: NO need to trim these reads as these are simulated reads
echo "# ------------------------------ #"
echo "XICRA join -i XICRA_analysis_PE --threads 4 --noTrim"
echo "..."
XICRA join -i XICRA_analysis_PE --threads 4 --noTrim | tee logs/XICRA_analysis_PE.join.log
echo ""

echo "# ------------------------------ #"
echo "XICRA miRNA -i XICRA_analysis_PE --threads 4 --software miraligner"
echo "..."
XICRA miRNA -i XICRA_analysis_PE --threads 4 --software miraligner | tee logs/XICRA_analysis_PE.miRNA.log
echo "" 

echo ""   
echo "`date` ... Paired End Test finished ..."
echo ""   
echo "# ------------------------------ #"
echo ""
echo ""
echo ""
echo ""
echo ""
echo ""
##########################################

##########################################
## test tRNA analysis
##########################################
echo ""   
echo "# ------------------------------ #"
echo "XICRA prep -i ./subset_tRNA/ -o XICRA_analysis_tRNA --single_end"
echo "..."
echo ""
XICRA prep -i ./subset_tRNA/ -o XICRA_analysis_tRNA --single_end | tee logs/XICRA_analysis_tRNA.prep.log
echo ""

echo "# ------------------------------ #"
echo "XICRA QC -i XICRA_analysis_tRNA --single_end --threads 4"
echo "..."
XICRA QC -i XICRA_analysis_tRNA --single_end --threads 4 | tee logs/XICRA_analysis_tRNA.qc.log
echo ""

# ATTENTION: NO need to trim these reads as these are simulated reads

echo "# ------------------------------ #"
echo "XICRA tRNA -i XICRA_analysis_tRNA --single_end --threads 4 --software mintmap"
echo "..."
XICRA tRNA -i XICRA_analysis_tRNA --noTrim --single_end --threads 4 --software mintmap | tee logs/XICRA_analysis_tRNA.tRNA.log
echo ""
echo "# ------------------------------ #"
echo ""
echo ""
echo ""   
echo "`date` ... tRNA Single End Test finished ..."


echo ""   
echo "`date` ... Finished ..."
##########################################
