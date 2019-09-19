## read data
# ---------------------------------------------------------------------------------
biotypes <- read.csv("/home/labs/lslab/jsanchez/DATA/XICRA/te", sep=",", stringsAsFactors=T, header=T, row.names = 1)
head(biotypes)

biotypes$biotype<- ifelse(rownames(biotypes)=="miRNA", "miRNA",
                          ifelse(rownames(biotypes)=="piRNA", "piRNA", 
                                 ifelse(rownames(biotypes)=="tRNA", "tRNA", 
                                        ifelse(rownames(biotypes)=="lincRNA", "lincRNA", 
                                              ifelse(rownames(biotypes)=="processed_transcript", "Processed Transcript", 
                                                 ifelse(rownames(biotypes)=="protein_coding", "Protein Coding", 
                                                        ifelse(rownames(biotypes)=="__alignment_not_unique", "Align NoUniq", 
                                                               ifelse(rownames(biotypes)=="__ambiguous", "Align Ambig", 
                                                                      ifelse(rownames(biotypes)=="__no_feature", "No Feature", 
                                                                             ifelse(rownames(biotypes)=="__not_aligned", "Not Align", 
                                                                                    ifelse(rownames(biotypes)=="antisense", "Antisense", 
                                                                                           ifelse(rownames(biotypes)=="misc_RNA", "misc_RNA", "Other"))))))))))))

## replace NA
biotypes[is.na(biotypes)] <- 0

# ---------------------------------------------------------------------------------
## summarize
sum.biotypes<- aggregate(biotypes[,-ncol(biotypes)], by=list(biotypes$biotype), "sum")
sum.biotypes$Group.1 <- factor(sum.biotypes$Group.1, levels=unique(biotypes$biotype))
View(sum.biotypes)

# ---------------------------------------------------------------------------------
## generate percentage
perc.biotypes<- cbind.data.frame(sum.biotypes$Group.1, 100 * prop.table( as.matrix(sum.biotypes[, 2:ncol(sum.biotypes)]), 2 ) )
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
