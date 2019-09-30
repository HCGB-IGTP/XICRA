## read data
# ---------------------------------------------------------------------------------
biotypes <- read.csv("/imppc/labs/lslab/share/data/proc_data/20190726_CTural_sRNAseq_ReAnalysis_170210/output_XICRA/counts_RNAbiotypes.csv", sep=",", stringsAsFactors=T, header=T, row.names = 1)
#biotypes <- read.csv("/imppc/labs/lslab/share/data/proc_data/20190301_MVives_Bioinfo_RNAseq_02/output_results/", sep=",", stringsAsFactors=T, header=T, row.names = 1)

head(biotypes)

## remove row
biotypes.remove <- c("total")
biotypes <- biotypes[!(row.names(biotypes) %in% biotypes.remove), ]


## types of RNA biotypes
rownames(biotypes)
#View(biotypes)

## get new names
biotypes$biotype<- ifelse(rownames(biotypes)=="miRNA", "miRNA",
                          ifelse(rownames(biotypes)=="piRNA", "piRNA", 
                                 ifelse(rownames(biotypes)=="tRNA", "tRNA", 
                                        ifelse(rownames(biotypes)=="lincRNA", "lincRNA", 
                                               ifelse(rownames(biotypes)=="rRNA", "rRNA", 
                                                  ifelse(rownames(biotypes)=="processed_transcript", "Processed Transcript", 
                                                     ifelse(rownames(biotypes)=="protein_coding", "Protein Coding", 
                                                          ifelse(rownames(biotypes)=="multimapping", "Align NoUniq", 
                                                                 ifelse(rownames(biotypes)=="Unassigned_Ambiguity", "No Feature", 
                                                                        ifelse(rownames(biotypes)=="unmapped", "Not Align", 
                                                                               ifelse(rownames(biotypes)=="antisense", "Antisense", 
                                                                                   ifelse(rownames(biotypes)=="misc_RNA", "misc_RNA", "Other"))))))))))))
library(ggplot2)
library(RColorBrewer)
library(reshape)
library(stringr)
## replace NA
biotypes[is.na(biotypes)] <- 0

# ---------------------------------------------------------------------------------
## summarize
sum.biotypes<- aggregate(biotypes[,-ncol(biotypes)], by=list(biotypes$biotype), "sum")
sum.biotypes$Group.1 <- factor(sum.biotypes$Group.1, levels=unique(biotypes$biotype))
#View(sum.biotypes)

# ---------------------------------------------------------------------------------
## generate percentage
perc.biotypes<- cbind.data.frame(sum.biotypes$Group.1, 100 * prop.table( as.matrix(sum.biotypes[, 2:ncol(sum.biotypes)]), 2 ) )
(perc.biotypes)

head(perc.biotypes)
names(perc.biotypes)[1]<- "Group.1"
# ---------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------
# Plot percentage Biotypes 
melt.biotypes <- melt(perc.biotypes)
melt.biotypes<- melt.biotypes[order(melt.biotypes$Group.1, melt.biotypes$value,  decreasing=T),]
ord<- rep(melt.biotypes$variable[1:144], 10)

head(melt.biotypes)
#ggplot(melt.biotypes, aes(x=reorder(variable, +value), y=value, fill=Group.1 )) +
#ggplot(melt.biotypes, aes(x=reorder(variable, ord), y=value, fill=Group.1 )) +

ggplot(melt.biotypes, aes(x=variable, y=value, fill=Group.1 )) +
  geom_bar(stat="identity") + scale_fill_brewer(palette="Spectral")


# ---------------------------------------------------------------------------------
# Plot PCA by percentage Biotypes -------------------------------------------------
perc.bio<- perc.biotypes
rownames(perc.bio)<- perc.bio$Group.1
perc.bio<- perc.bio[,-1]
perc.bio<- as.matrix(perc.bio)
perc.bio<- perc.bio[, order(perc.bio["Align NoUniq",])]
#perc.bio<- perc.bio[, order(perc.bio["miRNA",])]

#View(perc.bio)
list_legend <- c(
  "Align NoUniq",
  "Not Align",  
  "No Feature", 
  "Other",
  "Processed Transcript", 
  "Protein Coding", 
  "Antisense",
  "miRNA",  
  "rRNA",  
  "lincRNA",
  "misc_RNA"
  ##"piRNA",
)

perc.bio<- perc.bio[ list_legend, ]
barplot(perc.bio, col= rev(brewer.pal(11, "Spectral") ), xaxt="n", ylab = "%", main = "Biotypes distribution sample")
legend("bottomright",
       leg=list_legend,
       cex=0.6, 
       fill=rev(brewer.pal(11, "Spectral") ), 
       horiz=F, 
       ncol = 1)


# ---------------------------------------------------------------------------------
# Plot percentage Biotypes sorted by Subtype --------------------------------------

load(file="/imppc/labs/lslab/share/data/proc_data/CTural_sRNAseq_170210/data/sample.data.Rdata")
sample.data

rownames(sample.data) <- gsub(".count", "", sample.data$Sample_Code)

#perc.bio<- perc.bio[, rownames(sample.data)]

colors<-ifelse(sample.data$Subtype=="FO NP" 
               | sample.data$Subtype=="FO PROG" 
               | sample.data$Subtype=="F1 NP" 
               | sample.data$Subtype=="F1 PROG" 
               | sample.data$Subtype=="F2 NP" 
               | sample.data$Subtype=="F2 PROG", "yellow", 
               ifelse(sample.data$Subtype=="ALT ELEVAT", "royalblue1", 
                      ifelse(sample.data$Subtype=="HIPERPLASSIA NODULAR", "orange", 
                             ifelse(sample.data$Subtype=="ALT NORMAL", "green3", 
                                    ifelse(sample.data$Subtype=="VHC + monoinfectats", "maroon3", "mediumpurple1")))))


## 173 -->144
## VIH+VHC+               num yellow samples: 47
## ALT ELEVAT:            num royalblue1 samples: 21
## HIPERPLASSIA NODULAR:  num orange samples:11
## ALT NORMAL:            num green samples: 22
## VHC + monoinfectats:   num magenta samples: 22
## Healthy:               num mediumpurple1 samples: 21
#par(xpd=T)
#par(xpd=F)

table(sample.data$Subtype)
rownames(sample.data)

#samp.short<- str_sub(colnames(perc.bio), -3,-1)
perc.bio_here <- perc.bio[, rownames(sample.data)]
barplot(perc.bio_here, col= rev(brewer.pal(11, "Spectral") ),
        xaxt="n", 
        ylab = "%", main = "Biotypes distribution sample")

###
rect(0, -5, 56.46, 0, col="yellow", border = NULL)
rect(56.46, -5, 81.7, 0, col="royalblue1", border = NULL)
rect(81.7, -5, 94.91, 0, col="orange", border = NULL)
rect(94.91, -5, 121.34, 0, col="green3", border = NULL)
rect(121.34, -5, 147.77, 0, col="maroon3", border = NULL)
rect(147.77, -5, 173, 0, col="mediumpurple1", border = NULL)

legend(x=2, y=55,
       title = "Biotype",
       leg=list_legend,
       cex=0.6, 
       fill=rev(brewer.pal(11, "Spectral") ), 
       horiz=F, 
       ncol = 1)

legend(x="topright",
       title = "Sample group",
       leg=c("VIH+VHC+", "ALT ELEV", "HYPERPLASIA", "ALT NORM", "VHCpos Mono", "Healthy"), 
       cex=0.7, 
       fill=c("yellow", "royalblue1", "orange", "green3", "maroon3", "mediumpurple1"), 
       horiz=F, ncol = 2)
# ---------------------------------------------------------------------------------


