R.version #4.1.1 
rm(list = ls())
gc()

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#                            Core Microbiome                                 #                                   
#                                                                            #
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
Dir <- "/Users/simonmorvan/Desktop/UdeM/Projet/Ch2/16S_Data_Analysis/03_RData/"
setwd(Dir)

##### Packages #####
library(MicEco)

load("ps_bac_merg.RData")

Bac_Venn_ASV <- ps_venn(ps.bac.merg,
                        group = "Burning_intensity",
                        labels = list(col = "gray12", cex=1.2),
                        weight=F,
                        quantities = list(type=c("counts"),col = "gray12"),
                        col = "black",
                        fill= c("#9cccd8","#ffb53f","#f26419","#931125"),
                        adjust_labels=T)
Bac_Venn_ASV

Bac_Venn_relAb <- ps_venn(ps.bac.merg,
                        group = "Burning_intensity",
                        labels = list(col = "gray12", cex=1.5),
                        weight=T,
                        quantities = list(type=c("percent"),col = "gray12"),
                        col = "black",
                        fill= c("#9cccd8","#ffb53f","#f26419","#931125"),
                        adjust_labels=T)
Bac_Venn_relAb


Bac_Venn <- ps_venn(ps.bac.merg,
                        group = "Burning_intensity",
                        labels = list(col = "gray12", cex=1.1),
                        weight=F,
                        quantities = list(labels = c("522\n(<1%)",
                                                     "294\n(<1%)",
                                                     "505\n(<1%)",
                                                     "474\n(<1%)",
                                                     "25\n(<1%)",
                                                     "36\n(<1%)",
                                                     "17\n(<1%)",
                                                     "17\n(<1%)",
                                                     "14\n(<1%)",
                                                     "21\n(<1%)",
                                                     "70\n(<1%)",
                                                     "58\n(<1%)",
                                                     "93\n(<1%)",
                                                     "68\n(<1%)",
                                                     "1802\n(95%)")),
                        col = "black",
                        fill= c("#9cccd8","#ffb53f","#f26419","#931125"),
                        adjust_labels=T)
Bac_Venn


## Check # 

NullBurn <- prune_samples(ps.merg.tree@sam_data$Burning_intensity == "0",ps.merg.tree)
NullBurn<- prune_taxa(taxa_sums(NullBurn)>0, NullBurn)
# 469 ASVs

LowBurn <- prune_samples(ps.merg.tree@sam_data$Burning_intensity == "1",ps.merg.tree)
LowBurn<- prune_taxa(taxa_sums(LowBurn)>0, LowBurn)
# 399 ASVs

MedBurn <- prune_samples(ps.merg.tree@sam_data$Burning_intensity == "2",ps.merg.tree)
MedBurn<- prune_taxa(taxa_sums(MedBurn)>0, MedBurn)
# 431 ASVs

HiBurn <- prune_samples(ps.merg.tree@sam_data$Burning_intensity == "3",ps.merg.tree)
HiBurn<- prune_taxa(taxa_sums(HiBurn)>0, HiBurn)
# 365 ASVs


core <- Reduce(intersect, list(taxa_names(NullBurn),
                               taxa_names(LowBurn),
                               taxa_names(MedBurn),
                               taxa_names(HiBurn)))
# 222 ASVs

ps.core <-  prune_taxa(taxa_names(ps.merg.tree) %in% core,ps.merg.tree)
core_ab <- as.data.frame(colSums(otu_table(ps.core)))
tot_core_ab <- sum(core_ab)

ps_ab <- as.data.frame(colSums(otu_table(ps.merg.tree)))
tot_ps_ab <- sum(ps_ab)


rel_ab_core <- tot_core_ab/tot_ps_ab


burn <- Reduce(intersect, list(taxa_names(LowBurn),
                               taxa_names(MedBurn),
                               taxa_names(HiBurn)))
