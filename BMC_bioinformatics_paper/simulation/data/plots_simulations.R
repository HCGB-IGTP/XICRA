library(ggpubr)

## XICRA simulation
setwd("./")

#############################################################
## HS2500: 100 replicates x 100 miRNA percentage differences 0
#############################################################
perc0_df_simulations = read.table("./simulations_results_percDiff-0_XICRA.simulations.csv", sep = ',', header = 1, stringsAsFactors = TRUE)
head(perc0_df_simulations)

## Filter values: New generation and not Observed
perc0_new_creation_entries <- perc0_df_simulations[ with (perc0_df_simulations, grepl("New", variant)),]
perc0_NotObserved_entries <- perc0_df_simulations[ with (perc0_df_simulations, grepl("Not", variant)),]
perc0_df_simulations_final <- perc0_df_simulations[!rownames(perc0_df_simulations) %in% rownames(perc0_new_creation_entries),]
perc0_df_simulations_final <- perc0_df_simulations_final[!rownames(perc0_df_simulations_final) %in% rownames(perc0_NotObserved_entries),]

## make sure no Sensibility == 0
print (perc0_df_simulations_final[perc0_df_simulations_final$S == 0,])


## plot them all
perc0_S_plot <- ggboxplot(perc0_df_simulations_final, x="soft", y="S", color="type_read", 
                          xlab ="Software", ylab="Sensitivity", 
                          legend.title="Read type")

perc0_P_plot <- ggboxplot(perc0_df_simulations_final, x="soft", y="P", color="type_read", 
                          xlab ="Software", ylab="Precision", 
                          legend.title="Read type")

ggarrange(perc0_S_plot, perc0_P_plot, ncol=2, nrow=1, common.legend = TRUE, legend="bottom")
ggarrange(perc0_S_plot, perc0_P_plot, ncol=1, nrow=2, common.legend = TRUE, legend="top")

sum_perc0_new_creation_entries <- aggregate(perc0_new_creation_entries$obs, by=list(Category=perc0_new_creation_entries$soft, perc0_new_creation_entries$type_read), FUN=sum)
sum_perc0_new_creation_entries
sum_perc0_new_bar <- ggbarplot(sum_perc0_new_creation_entries, x="Category", y="x", fill="Group.2", 
                               xlab ="Software", ylab="New generation", position = position_dodge(0.9),
                               legend.title="Read type")

ggarrange(perc0_S_plot + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                               axis.title.x = element_blank()) + rremove('x.axis'),
          perc0_P_plot + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                               axis.title.x = element_blank()) + rremove('x.axis'),
          sum_perc0_new_bar,
          ncol=1, nrow=3, common.legend = TRUE, legend="top", align = "v")



#############################################################

#############################################################
## HS2500: 100 replicates x 100 miRNA percentage differences 8
#############################################################
perc8_df_simulations = read.table("./simulations_results_percDiff-8_XICRA.simulations.csv", sep = ',', header = 1, stringsAsFactors = TRUE)
head(perc8_df_simulations)

## Filter values: New generation and not Observed
perc8_new_creation_entries <- perc8_df_simulations[ with (perc8_df_simulations, grepl("New", variant)),]
perc8_NotObserved_entries <- perc8_df_simulations[ with (perc8_df_simulations, grepl("Not", variant)),]
perc8_df_simulations_final <- perc8_df_simulations[!rownames(perc8_df_simulations) %in% rownames(perc8_new_creation_entries),]
perc8_df_simulations_final <- perc8_df_simulations_final[!rownames(perc8_df_simulations_final) %in% rownames(perc8_NotObserved_entries),]

## make sure no Sensibility == 0
print (perc8_df_simulations_final[perc8_df_simulations_final$S == 0,])


## plot them all
perc8_S_plot <- ggboxplot(perc8_df_simulations_final, x="soft", y="S", color="type_read", 
                         xlab ="Software", ylab="Sensitivity", 
                         legend.title="Read type")

perc8_P_plot <- ggboxplot(perc8_df_simulations_final, x="soft", y="P", color="type_read", 
                         xlab ="Software", ylab="Precision", 
                         legend.title="Read type")

ggarrange(perc8_S_plot, perc8_P_plot, ncol=2, nrow=1, common.legend = TRUE, legend="bottom")
ggarrange(perc8_S_plot, perc8_P_plot, ncol=1, nrow=2, common.legend = TRUE, legend="top")

sum_perc8_new_creation_entries <- aggregate(perc8_new_creation_entries$obs, 
                                            by=list(Category=perc8_new_creation_entries$soft, 
                                                    perc8_new_creation_entries$type_read), FUN=sum)
sum_perc8_new_creation_entries
sum_perc8_new_bar <- ggbarplot(sum_perc8_new_creation_entries, x="Category", y="x", fill="Group.2", 
                               xlab ="Software", ylab="New generation", position = position_dodge(0.9),
                               legend.title="Read type")

ggarrange(perc8_S_plot + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                               axis.title.x = element_blank()) + rremove('x.axis'),
          perc8_P_plot + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                               axis.title.x = element_blank()) + rremove('x.axis'),
          sum_perc8_new_bar,
          ncol=1, nrow=3, common.legend = TRUE, legend="top", align = "v")

#############################################################

#############################################################
##### MiSeq V1 profile: 10 replicates x 100 miRNA
#############################################################
MSv1_df_simulations = read.table("./simulations_MSv1.csv", sep = ',', header = 1, stringsAsFactors = TRUE)
head(MSv1_df_simulations)

## Filter values: New generation and not Observed
MSv1_new_creation_entries <- MSv1_df_simulations[ with (MSv1_df_simulations, grepl("New", variant)),]
MSv1_NotObserved_entries <- MSv1_df_simulations[ with (MSv1_df_simulations, grepl("Not", variant)),]
MSv1_df_simulations_final <- MSv1_df_simulations[!rownames(MSv1_df_simulations) %in% rownames(MSv1_new_creation_entries),]
MSv1_df_simulations_final <- MSv1_df_simulations_final[!rownames(MSv1_df_simulations_final) %in% rownames(MSv1_NotObserved_entries),]

## make sure no Sensibility == 0
print (MSv1_df_simulations_final[MSv1_df_simulations_final$S == 0,])


## plot them all
MSv1_S_plot <- ggboxplot(MSv1_df_simulations_final, x="soft", y="S", color="type_read", 
          xlab ="Software", ylab="Sensitivity", 
          legend.title="Read type")

MSv1_P_plot <- ggboxplot(MSv1_df_simulations_final, x="soft", y="P", color="type_read", 
          xlab ="Software", ylab="Precision", 
          legend.title="Read type")

ggarrange(MSv1_S_plot, MSv1_P_plot, ncol=2, nrow=1, common.legend = TRUE, legend="bottom")
ggarrange(MSv1_S_plot, MSv1_P_plot, ncol=1, nrow=2, common.legend = TRUE, legend="top")

sum_MSv1_new_creation_entries <- aggregate(MSv1_new_creation_entries$obs, by=list(Category=MSv1_new_creation_entries$soft, MSv1_new_creation_entries$type_read), FUN=sum)
sum_MSv1_new_bar <- ggbarplot(sum_MSv1_new_creation_entries, x="Category", y="x", fill="Group.2", 
                              xlab ="Software", ylab="New generation", position = position_dodge(0.9),
                              legend.title="Read type")


#############################################################
