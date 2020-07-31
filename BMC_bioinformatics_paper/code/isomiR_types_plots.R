devtools::install_github("HCGB-IGTP/XICRA.stats")
tmp_folder <- "/home/jfsanchez/test_gmt"
setwd("/home/jfsanchez/git_repos/XICRA/BMC_bioinformatics_paper/")

library(XICRA.stats)
library(reshape2)
library(ggplot2)

## save preRank files for gsea
save_preRank_file <- function(data_df, path2storeResults, tag){
  # save preranked input files for input into GSEA jave applet tool
  write.table(data_df[,c("UID","log2FoldChange")], 
              file=file.path(path2storeResults, paste0(tag, ".data.rnk")), 
              sep="\t", row.names=F, col.names=F, quote=F)
}

## get df items from colum
get_items <- function(list, df_frame, col_name) {
  v <- c()
  for (p in list) {
    item <- df_frame[p,][[col_name]]
    print (item)
    v <- c(v, item)
  }
  return(v)
}

## calculate stats for simulations
calculates_Stats <- function(softname, data_df_percDiff0, data_df_percDiff8, all_isomir_types=FALSE) {
  ## perc_diff0
  data_df_percDiff8 <- data_df_percDiff8[data_df_percDiff8$soft==softname,]
  data_df_percDiff8$variant <- as.character(data_df_percDiff8$variant)
  data_df_percDiff8$variant[data_df_percDiff8$variant == "NA"] <- "Canonical"
  data_df_percDiff8$variant[is.na(data_df_percDiff8$variant)] <- "Canonical"  
  data_df_percDiff8 <- tidyr::separate(data_df_percDiff8, variant, c('type', 'int'), sep = ':', remove=FALSE)
  data_df_percDiff8$type <- as.character(data_df_percDiff8$type)
  data_df_given_PE_percDiff8 <- data_df_percDiff8[data_df_percDiff8$type_read=='PE',]
  
  ## calculate geneSet stats
  PE_8_stats <- XICRA.stats::create_geneSet_stats_simulations(data_df_given_PE_percDiff8)
  
  ## perc_diff0
  data_df_percDiff0 <- data_df_percDiff0[data_df_percDiff0$soft==softname,]
  data_df_percDiff0$variant <- as.character(data_df_percDiff0$variant)
  data_df_percDiff0$variant[data_df_percDiff0$variant == "NA"] <- "Canonical"
  data_df_percDiff0$variant[is.na(data_df_percDiff0$variant)] <- "Canonical"  
  data_df_percDiff0 <- tidyr::separate(data_df_percDiff0, variant, c('type', 'int'), sep = ':', remove=FALSE)
  data_df_percDiff0$type <- as.character(data_df_percDiff0$type)
  
  data_df_given_PE_percDiff0 <- data_df_percDiff0[data_df_percDiff0$type_read=='PE',]
  data_df_given_R1_percDiff0 <- data_df_percDiff0[data_df_percDiff0$type_read=='R1',]
  data_df_given_R2_percDiff0 <- data_df_percDiff0[data_df_percDiff0$type_read=='R2',]
  
  ####
  PE_0_stats <- XICRA.stats::create_geneSet_stats_simulations(data_df_given_PE_percDiff0)
  R1_stats <- XICRA.stats::create_geneSet_stats_simulations(data_df_given_R1_percDiff0)
  R2_stats <- XICRA.stats::create_geneSet_stats_simulations(data_df_given_R2_percDiff0)
  
  #isotype plots per method
  counts.longData <- calculate_average_readcounts_plot(PE_0_df = PE_0_stats, PE_8_df = PE_8_stats,
                                    SE_R1_df = R1_stats, SE_R2_df = R2_stats, use_all = all_isomir_types)
  
  UID_longData <- calculate_unique_readCounts_plot(PE_0_df = PE_0_stats, PE_8_df = PE_8_stats,
                                                       SE_R1_df = R1_stats, SE_R2_df = R2_stats, use_all = all_isomir_types)
  
  #return
  list2return <- list("counts.longData" =  counts.longData,
                      "UID_longData" = UID_longData)
  return(list2return)
}

## plot average
average_read_plot <- function(counts.longData) {
  p <- (ggplot(counts.longData, aes(x = variable, y = iso_types)) + 
          geom_raster(aes(fill=value)) + 
          scale_fill_gradient(low="grey90", high="red") +
          labs(x="Sequencing mode", y="Isomir types", title="Average counts per isomir type") +
          theme_bw() +
          geom_point(aes(size=value)) +
          theme_bw() + 
          theme(axis.text.x=element_text(size=12, angle=90, vjust=0.3),
                axis.text.y=element_text(size=12),
                plot.title=element_text(size=14),
                legend.text=element_text(size=12)))
  return(p)
}

## plot unique 
unique_UID_plot <-function(UID_longData) {
  print(ggplot(UID_longData, aes(x = variable, y = iso_types)) + 
          geom_raster(aes(fill=value)) + 
          scale_fill_gradient(low="grey90", high="red") +
          labs(x="Sequencing mode", y="Isomir types", title="Unique isomiRs") +
          theme_bw() +
          geom_point(aes(size=value)) +
          theme_bw() + 
          theme(axis.text.x=element_text(size=12, angle=90, vjust=0.3),
                axis.text.y=element_text(size=12),
                plot.title=element_text(size=14),
                legend.text=element_text(size=12)))
  
}

## calculate data for average_read counts
calculate_average_readcounts_plot <- function (PE_0_df, PE_8_df, SE_R1_df, SE_R2_df, use_all=TRUE) {
  if (use_all) {
    ## or use them all
    union1 <- union(SE_R1_df$class, SE_R2_df$class)
    union2 <- union(PE_0_df$class, union1)
    union3 <- union(PE_8_df$class, union2)
    length(union3)
    
  } else {
    ## use 10 items
    union3 <- PE_0_df$class[c(1:10)]
  }
  
  stats.table.counts<-data.frame("PE_0"=c(1:length(union3)), "PE_8"=c(1:length(union3)), "SR1"=c(1:length(union3)),"SR2"=c(1:length(union3)))
  stats.table.counts$iso_types<-union3
  
  stats.table.counts[,1]<-get_items(stats.table.counts$iso_types, PE_0_df, "average.read.counts")
  stats.table.counts[,2]<-get_items(stats.table.counts$iso_types, PE_8_df, "average.read.counts")
  stats.table.counts[,3]<-get_items(stats.table.counts$iso_types, SE_R1_df, "average.read.counts")
  stats.table.counts[,4]<-get_items(stats.table.counts$iso_types, SE_R2_df, "average.read.counts")
  stats.table.counts[is.na(stats.table.counts)] <- 0
  
  counts.longData<-melt(stats.table.counts)
  counts.longData<-counts.longData[counts.longData$value!=0,]
  
  ## plot
  average_read_plot(counts.longData)
  
  return(counts.longData)
}

## calculate data for unique read counts plot
calculate_unique_readCounts_plot <- function(PE_0_df, PE_8_df, SE_R1_df, SE_R2_df, use_all=TRUE) {
  if (use_all) {
    ## or use them all
    union1 <- union(SE_R1_df$class, SE_R2_df$class)
    union2 <- union(PE_0_df$class, union1)
    union3 <- union(PE_8_df$class, union2)
    length(union3)
    
  } else {
    ## use 10 items
    union3 <- PE_0_df$class[c(1:10)]
    
  }
  #unique ID plots per method
  stats.table_UID<-data.frame("PE_0"=c(1:length(union3)), "PE_8"=c(1:length(union3)), "SR1"=c(1:length(union3)),"SR2"=c(1:length(union3)))
  stats.table_UID$iso_types<-union3
  
  stats.table_UID[,1]<- get_items(stats.table_UID$iso_types, PE_0_df, "n.unique")
  stats.table_UID[,2]<- get_items(stats.table_UID$iso_types, PE_8_df, "n.unique")
  stats.table_UID[,3]<- get_items(stats.table_UID$iso_types, SE_R1_df, "n.unique")
  stats.table_UID[,4]<- get_items(stats.table_UID$iso_types, SE_R2_df, "n.unique")
  stats.table_UID[is.na(stats.table_UID)] <- 0
  
  UID_longData<-melt(stats.table_UID)
  UID_longData<-UID_longData[UID_longData$value!=0,]
  
  ## plot
  unique_UID_plot(UID_longData)
  
  return(UID_longData)
}


#######################################
### BST SAMPLES
#######################################

## files
PE_0 <- file.path(tmp_folder, "PE_0_DESeq2_table.tsv")
PE_8 <- file.path(tmp_folder, "PE_8_DESeq2_table.tsv")
SE_R1 <- file.path(tmp_folder, "SE_R1_DESeq2_table.tsv")
SE_R2 <- file.path(tmp_folder, "SE_R2_DESeq2_table.tsv")

## read data
PE_0_data <- XICRA.stats::prepare_data(PE_0, "DESeq2")
PE_8_data <- XICRA.stats::prepare_data(PE_8, "DESeq2")
SE_R1_data <- XICRA.stats::prepare_data(SE_R1, "DESeq2")
SE_R2_data <- XICRA.stats::prepare_data(SE_R2, "DESeq2")

## Prepare genet set data
PE_0_data_geneSet <- XICRA.stats::build_geneSet_collection_DESeq2(PE_0_data, "PE_0", "./")
PE_8_data_geneSet <- XICRA.stats::build_geneSet_collection_DESeq2(PE_8_data, "PE_8", "./")
SE_R1_data_geneSet <- XICRA.stats::build_geneSet_collection_DESeq2(SE_R1_data, "SE_R1", "./")
SE_R2_data_geneSet <- XICRA.stats::build_geneSet_collection_DESeq2(SE_R2_data, "SE_R2", "./")

##
save_preRank_file(data_df = PE_0_data, path2storeResults = "./", tag = "PE_0")
save_preRank_file(data_df = PE_8_data, path2storeResults = "./", tag = "PE_8")
save_preRank_file(data_df = SE_R1_data, path2storeResults = "./", tag = "SE_R1")
save_preRank_file(data_df = SE_R2_data, path2storeResults = "./", tag = "SE_R2")

## do gsea analysis

#isotype plots per method
#isotype plots per method
calculate_average_readcounts_plot(PE_0_df = PE_0_data_geneSet$geneSet_stats, PE_8_df = PE_8_data_geneSet$geneSet_stats,
                                  SE_R1_df = SE_R1_data_geneSet$geneSet_stats, SE_R2_df = SE_R2_data_geneSet$geneSet_stats,
                                  use_all = TRUE)

calculate_unique_readCounts_plot(PE_0_df = PE_0_data_geneSet$geneSet_stats, PE_8_df = PE_8_data_geneSet$geneSet_stats,
                                 SE_R1_df = SE_R1_data_geneSet$geneSet_stats, SE_R2_df = SE_R2_data_geneSet$geneSet_stats,
                                 use_all = TRUE)

#######################################

#######################################
## GSE114923
#######################################
## files
PE_0_GSE114923 <- "./analysis_GSE114923/results/PE_0_DESeq2_table.tsv"
PE_8_GSE114923 <- "./analysis_GSE114923/results/PE_8_DESeq2_table.tsv"
SE_R1_GSE114923 <- "./analysis_GSE114923/results/SE_R1_DESeq2_table.tsv"
SE_R2_GSE114923 <- "./analysis_GSE114923/results/SE_R2_DESeq2_table.tsv"

## read data
PE_0_GSE114923_data <- XICRA.stats::prepare_data(PE_0_GSE114923, "DESeq2")
PE_8_GSE114923_data <- XICRA.stats::prepare_data(PE_8_GSE114923, "DESeq2")
SE_R1_GSE114923_data <- XICRA.stats::prepare_data(SE_R1_GSE114923, "DESeq2")
SE_R2_GSE114923_data <- XICRA.stats::prepare_data(SE_R2_GSE114923, "DESeq2")

## Prepare genet set data
PE_0_data_geneSet_GSE114923 <- XICRA.stats::build_geneSet_collection_DESeq2(PE_0_GSE114923_data, "PE_0", tmp_folder)
PE_8_data_geneSet_GSE114923 <- XICRA.stats::build_geneSet_collection_DESeq2(PE_8_GSE114923_data, "PE_8", tmp_folder)
SE_R1_data_geneSet_GSE114923 <- XICRA.stats::build_geneSet_collection_DESeq2(SE_R1_GSE114923_data, "SE_R1", tmp_folder)
SE_R2_data_geneSet_GSE114923 <- XICRA.stats::build_geneSet_collection_DESeq2(SE_R2_GSE114923_data, "SE_R2", tmp_folder)

##
save_preRank_file(data_df = PE_0_data, path2storeResults = "./", tag = "PE_0")
save_preRank_file(data_df = PE_8_data, path2storeResults = "./", tag = "PE_8")
save_preRank_file(data_df = SE_R1_data, path2storeResults = "./", tag = "SE_R1")
save_preRank_file(data_df = SE_R2_data, path2storeResults = "./", tag = "SE_R2")

## do gsea analysis

#isotype plots per method
calculate_average_readcounts_plot(PE_0_df = PE_0_data_geneSet_GSE114923$geneSet_stats, PE_8_df = PE_8_data_geneSet_GSE114923$geneSet_stats,
                                  SE_R1_df = SE_R1_data_geneSet_GSE114923$geneSet_stats, SE_R2_df = SE_R2_data_geneSet_GSE114923$geneSet_stats,
                                  use_all = TRUE)

calculate_unique_readCounts_plot(PE_0_df = PE_0_data_geneSet_GSE114923$geneSet_stats, PE_8_df = PE_8_data_geneSet_GSE114923$geneSet_stats,
                                  SE_R1_df = SE_R1_data_geneSet_GSE114923$geneSet_stats, SE_R2_df = SE_R2_data_geneSet_GSE114923$geneSet_stats,
                                  use_all = TRUE)
#######################################


#######################################
## Simulations
#######################################

## read simulations
data_df_given_percDiff8 <- read.table("./simulation/data/simulations_results_percDiff-8_XICRA.simulations.csv", sep=",", header=1)
data_df_given_percDiff8 <- data_df_given_percDiff8[!grepl("NotOb", data_df_given_percDiff8$variant),]

data_df_given_percDiff0 <- read.table("./simulation/data/simulations_results_percDiff-0_XICRA.simulations.csv", sep=",", header=1)
data_df_given_percDiff0 <- data_df_given_percDiff0[!grepl("NotOb", data_df_given_percDiff0$variant),]

## miraligner
miraligner_data <- calculates_Stats("miraligner", data_df_given_percDiff0, data_df_given_percDiff8, all_isomir_types = TRUE)
head(miraligner_data$counts.longData)
head(miraligner_data$UID_longData)

average_read_plot(miraligner_data$counts.longData)
unique_UID_plot(miraligner_data$UID_longData)

## sRNAbench
sRNAbench_data <- calculates_Stats("sRNAbench", data_df_given_percDiff0, data_df_given_percDiff8, all_isomir_types = TRUE)
head(sRNAbench_data$counts.longData)
head(sRNAbench_data$UID_longData)

average_read_plot(sRNAbench_data$counts.longData)
unique_UID_plot(sRNAbench_data$UID_longData)

## optimir
optimir_data <- calculates_Stats("optimir", data_df_given_percDiff0, data_df_given_percDiff8, all_isomir_types = TRUE)
head(optimir_data$counts.longData)
head(optimir_data$UID_longData)

average_read_plot(optimir_data$counts.longData)
unique_UID_plot(optimir_data$UID_longData)
