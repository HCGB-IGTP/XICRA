## install XICRA.stats
devtools::install_github("HCGB-IGTP/XICRA.stats")

## install HCGB.IGTP.Danalysis
devtools::install_github("HCGB-IGTP/HCGB.IGTP.Danalysis")

## load modules
library(XICRA.stats)
library(HCGB.IGTP.DAnalysis)
library(DESeq2)
library("RColorBrewer")
library(ggplot2)
library(gtools)
library(pheatmap)
library("BiocParallel")
library("apeglm")
library(vsn)
library(hexbin)
library("gplots")

## exploratory
exploratory_plots <- function(dds_object, OUTPUT_dir){
  ############################
  # Exploratory plots 
  ############################
  
  # Dispersion plot 
  jpeg(file.path(OUTPUT_dir, "general_dispersion_plot.jpeg"))
  plotDispEsts(dds_object)
  dev.off()
  
  # Top 50 genes:
  select <- order(rowMeans(counts(dds_object,normalized=TRUE)),decreasing=TRUE)[1:50]
  hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(1000)
  vsd <- varianceStabilizingTransformation(dds_object, blind = FALSE)
  
  ##order by condition
  df_condition <- as.data.frame(colData(dds_object)["condition"])
  colorsVector_condition = ifelse(df_condition["condition"]=="case","red", "blue")
  sorted_condition <- mixedorder(dds_object$condition)
  pdf(file.path(OUTPUT_dir, "general_heatmap_top50_byCondition.pdf"),
      pointsize = 8)
  heatmap.2(assay(vsd)[select,sorted_condition], 
            col = hmcol,
            Rowv = FALSE, Colv = FALSE, scale="none", 
            ColSideColors=colorsVector_condition[sorted_condition],
            dendrogram="none", trace="none", margin=c(15, 12))
  par(lend = 1)           # square line ends for the color legend
  legend("topright", legend=c("Case", "Control"),
         col=c("red", "blue"), 
         lty=1, title="Condition")
  dev.off()
  
  
  # Clustering:
  # Get sample-to-sample distances
  distsRL <- dist(t(assay(vsd)))
  mat <- as.matrix(distsRL)
  order_matrix <- match(colnames(mat), dds_object$sample_name)
  pdf(file.path(OUTPUT_dir,"general_heatmap_samplesDistance.pdf"))
  heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(15, 12), 
            ColSideColors=colorsVector_condition)
  legend("topright", legend=c("Case", "Control"),
         col=c("red", "blue"), 
         lty=1, title="Condition")
  dev.off()
  
  ## PCA
  pdf(file.path(OUTPUT_dir,"PCA_condition.pdf"))
  plt <- plotPCA(vsd, intgroup="condition")
  print(plt)
  dev.off()
  
  ### cooks distance
  pdf(file.path(OUTPUT_dir,"cooks_distance.pdf"))
  par(mar=c(8,5,2,2))
  boxplot(log10(assays(dds_object)[["cooks"]]), range=0, las=2)
  dev.off()
  
}

## DE analysis condition
analysis_DEseq <- function(counts_file, target_file, outdir, int_threads) {
  ## Create list object 
  object_DESeq <- list(
    "counts"=HCGB.IGTP.DAnalysis::discard_0_counts(counts = counts_file), 
    "target"=target_file)
  
  ## adjust samples in data and in sample sheet
  object_DESeq$counts <- object_DESeq$counts[, colnames(object_DESeq$counts) %in% rownames(object_DESeq$target) ]
  object_DESeq$counts <- object_DESeq$counts[, sort(colnames(object_DESeq$counts)) ]
  
  ##
  logical_vec <- rownames(object_DESeq$target) %in% colnames(object_DESeq$counts)
  object_DESeq$target <- object_DESeq$target[logical_vec,]
  match(rownames(object_DESeq$target), colnames(object_DESeq$counts)) 
  
  OUTPUT_Data_dir_test <- file.path(OUTPUT_Data_dir, outdir)
  dir.create(OUTPUT_Data_dir_test)
  
  ## comparison1 Condition: Chron vs. Control
  ddsFullCountTable_Condition <- DESeq2::DESeqDataSetFromMatrix( 
    countData = object_DESeq$counts, 
    colData = object_DESeq$target, 
    design = formula(~condition) )
  
  ddsFullCountTable_Condition$condition<- relevel(ddsFullCountTable_Condition$condition, ref="control")
  dds_Condition <- DESeq(ddsFullCountTable_Condition, parallel = TRUE)
  
  ## exploratory dds_Condition
  exploratory_plots_dir <- file.path(OUTPUT_Data_dir_test, "exploratory_condition")
  dir.create(exploratory_plots_dir)
  exploratory_plots(dds_object = dds_Condition, OUTPUT_dir = exploratory_plots_dir)
  
  ## check 
  target_file['samples'] <- NULL
  print("resultsNames(dds_Condition)") 
  print(resultsNames(dds_Condition))
  res_Case_vs_Contrl = try(HCGB.IGTP.DAnalysis::DESeq2_HCGB_function(dds_object = dds_Condition, coef_n = 2, 
                                                                     name= "comparison1_Condition", 
                                                                     numerator = "case", denominator = "control", 
                                                                     OUTPUT_Data_dir = OUTPUT_Data_dir_test, 
                                                                     df_treatment_Ind = target_file,
                                                                     threads = as.numeric(int_threads)))
}

###############################
## start analysis
###############################

## rscript
setwd("/imppc/labs/lslab/share/data/proc_data/20200714_XICRA_analysis_BMC_Bioinformatics/data/NCBI_projects/final_data_BMC/DE_analysis")

OUTPUT_Data_dir <- "./results/"
dir.create(OUTPUT_Data_dir)

##################
## sample sheet
##################
## simple sample sheet
samples <- c('SRR7218579', 'SRR7218580', 'SRR7218581', 'SRR7218582', 'SRR7218583', 'SRR7218584', 'SRR7218585', 'SRR7218586')
condition <- c('control', 'case', 'case', 'case', 'case', 'control', 'control', 'control')
de_sample_sheet <-  data.frame(row.names=samples, condition=condition, samples=samples)

##################
## PE data
##################
## PE miraligner
data_PE_miraligner <- "./data/PE/miRNA_expression-miraligner.csv"
data_PE_miraligner_data <- XICRA.stats::parse_XICRA(data_PE_miraligner)

analysis_DEseq(data_PE_miraligner_data$isomir_data, de_sample_sheet, outdir = 'PE_0_isomir-miraligner', int_threads = 2)
analysis_DEseq(data_PE_miraligner_data$miRNA_data, de_sample_sheet, outdir = 'PE_0_mirna-miraligner', int_threads = 2)

##################
## PE data
##################
## PE miraligner
data_PE_8_miraligner <- "./data/PE_8/miRNA_expression-miraligner.csv"
data_PE_8_miraligner_data <- XICRA.stats::parse_XICRA(data_PE_8_miraligner)

analysis_DEseq(data_PE_8_miraligner_data$isomir_data, de_sample_sheet, outdir = 'PE_8_isomir-miraligner', int_threads = 2)
analysis_DEseq(data_PE_8_miraligner_data$miRNA_data, de_sample_sheet, outdir = 'PE_8_mirna-miraligner', int_threads = 2)


##################
## SE R1 data
##################
## SE miraligner
data_SE_miraligner <- "./data/SE/miRNA_expression-miraligner.csv"
data_SE_miraligner_data <- XICRA.stats::parse_XICRA(data_SE_miraligner)

analysis_DEseq(data_SE_miraligner_data$isomir_data, de_sample_sheet, outdir = 'SE_isomir-miraligner', int_threads = 2)
analysis_DEseq(data_SE_miraligner_data$miRNA_data, de_sample_sheet, outdir = 'SE_mirna-miraligner', int_threads = 2)

##################
## SE R2 data
##################
## SE miraligner
data_SE_R2_miraligner <- "./data/SE_R2/miRNA_expression-miraligner.csv"
data_SE_R2_miraligner_data <- XICRA.stats::parse_XICRA(data_SE_R2_miraligner)

analysis_DEseq(data_SE_R2_miraligner_data$isomir_data, de_sample_sheet, outdir = 'SE_R2_isomir-miraligner', int_threads = 2)
analysis_DEseq(data_SE_R2_miraligner_data$miRNA_data, de_sample_sheet, outdir = 'SE_R2_mirna-miraligner', int_threads = 2)


##################################
## sessionInfo details
##################################
sessionInfo()

