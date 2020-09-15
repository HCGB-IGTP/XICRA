#!/usr/bin/env Rscript
library("optparse")

## get options
option_list = list(
  make_option(c("-f", "--file"), type="character", help="biotypes file", metavar="character"),
  make_option(c("-o", "--output"), type="character", help="name", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

## get message
if (is.null(opt$file)){
  print_help(opt_parser)
  stop("No arguments provided", call.=FALSE)
}

## read data
# ---------------------------------------------------------------------------------
## miRNAseq
biotypes <- read.csv(opt$file, sep=",", stringsAsFactors=T, header=T, row.names = 1)

## remove row
biotypes.remove <- c("total")
biotypes <- biotypes[!(row.names(biotypes) %in% biotypes.remove), ]

## types of RNA biotypes
#rownames(biotypes)

## get new names
biotypes$biotype<- ifelse(rownames(biotypes)=="miRNA", "miRNA",
                          ifelse(rownames(biotypes)=="piRNA", "piRNA", 
                                 ifelse(rownames(biotypes)=="tRNA", "tRNA", 
                                        ifelse(rownames(biotypes)=="lincRNA", "lincRNA", 
                                               ifelse(rownames(biotypes)=="rRNA", "rRNA", 
                                                  ifelse(rownames(biotypes)=="processed_transcript", "Processed Transcript", 
                                                     ifelse(rownames(biotypes)=="protein_coding", "Protein Coding", 
                                                          ifelse(rownames(biotypes)=="Unassigned_MultiMapping", "Align NoUniq", 
                                                                 ifelse(rownames(biotypes)=="multimapping", "Align NoUniq", 
                                                                      ifelse(rownames(biotypes)=="Unassigned_Ambiguity", "No Feature", 
                                                                        ifelse(rownames(biotypes)=="Unassigned_NoFeatures", "No Feature", 
                                                                          ifelse(rownames(biotypes)=="unmapped", "Not Align", 
                                                                                 ifelse(rownames(biotypes)=="antisense", "Antisense", 
                                                                                     ifelse(rownames(biotypes)=="misc_RNA", "misc_RNA", "Other"))))))))))))))
library(ggplot2)
library(RColorBrewer)
library(reshape)
library(stringr)

## replace NA
biotypes[is.na(biotypes)] <- 0

# ---------------------------------------------------------------------------------
## summarize
# ---------------------------------------------------------------------------------
sum.biotypes<- aggregate(biotypes[,-ncol(biotypes)], by=list(biotypes$biotype), "sum")
sum.biotypes$Group.1 <- factor(sum.biotypes$Group.1, levels=unique(biotypes$biotype))

## generate percentage
perc.biotypes<- cbind.data.frame(sum.biotypes$Group.1, 100 * prop.table( as.matrix(sum.biotypes[, 2:ncol(sum.biotypes)]), 2 ) )
names(perc.biotypes)[1]<- "Biotypes"

# ---------------------------------------------------------------------------------
# Plot percentage Biotypes ggplot
# ---------------------------------------------------------------------------------
melt.biotypes <- melt(perc.biotypes)
melt.biotypes<- melt.biotypes[order(melt.biotypes$Biotypes, melt.biotypes$value,  decreasing=T),]
names(melt.biotypes)[2] <- 'Samples'
names(melt.biotypes)[3] <- 'Value'

p <- ggplot(melt.biotypes, aes(x=Samples, y=Value, fill=Biotypes)) +
  geom_bar(stat="identity", colour="black") + 
  scale_fill_brewer(palette="Set3") + theme_classic() + 
  theme(axis.text.x = element_text(face = "bold", size = 10, angle = 90))
  
## save in pdf
ggsave(opt$output, p, width = 11, height = 8.5)

  
# ---------------------------------------------------------------------------------
# Plot percentage Biotypes basic R
# ---------------------------------------------------------------------------------
# perc.bio<- perc.biotypes
# rownames(perc.bio)<- perc.bio$Group.1
# perc.bio<- perc.bio[,-1]
# perc.bio<- as.matrix(perc.bio)
# perc.bio<- perc.bio[, order(perc.bio["Align NoUniq",])]
# #perc.bio<- perc.bio[, order(perc.bio["miRNA",])]
# 
# #View(perc.bio)
# list_legend <- c(
#   #"Align NoUniq",
#   "Not Align",  
#   "No Feature", 
#   "Other",
#   "Processed Transcript", 
#   "Protein Coding", 
#   "Antisense",
#   "miRNA",  
#   "rRNA",  
#   "lincRNA",
#   "misc_RNA"
#   ##"piRNA",
# )
# 
# perc.bio<- perc.bio[ list_legend, ]
# barplot_biotypes <- barplot(perc.bio, col= rev(brewer.pal(11, "Spectral") ), 
#                             ylab = "%", main = "Biotypes distribution sample")
# legend("bottomright",
#        leg=list_legend,
#        cex=0.6, inset=c(-0.15,0),
#        xpd=TRUE,
#        fill=rev(brewer.pal(11, "Spectral") ), 
#        horiz=F, 
#        ncol = 1)

