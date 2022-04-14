#### GETTING READY ####
R.version #4.1.1 
rm(list = ls())
gc()

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#                       Indicative species  script                           #                                   
#                                                                            #                   
#                                                                            #
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#


setwd("/Users/simonmorvan/Desktop/UdeM/Projet/Ch2/ITS_Data_Analysis/03_RData/")

##### Packages #####
library(phyloseq)
library(vegan)
library(indicspecies)

##### Data #####
load("ps.merg.tree.RData")

ASV_ab <- as.data.frame(otu_table(ps.merg.tree))
Tax_tab <- as.data.frame(tax_table(ps.merg.tree))
Sdata <- as.data.frame(sample_data(ps.merg.tree))


ASV_fun_indval <- multipatt(ASV_ab,Sdata$Burning_intensity,
                               control=how(nperm=999), 
                               duleg =TRUE,
                               print.perm =TRUE)

summary(ASV_fun_indval)
P_adj <- ASV_fun_indval$sign
P_adj$p_val_adj_holm <- p.adjust(P_adj$p.value, "holm")
P_adj$p_val_adj_fdr <- p.adjust(P_adj$p.value, "fdr")


# No significant indicative species 

##### MULTIPATT comb sites ####

ASV_fun_indval_comb <- multipatt(ASV_ab,Sdata$Burning_intensity,
                            control=how(nperm=999), 
                            duleg =FALSE,
                            print.perm =TRUE)

P_adj_comb <- ASV_fun_indval_comb$sign
P_adj_comb$p_val_adj_holm <- p.adjust(P_adj_comb$p.value, "holm")
P_adj_comb$p_val_adj_fdr <- p.adjust(P_adj_comb$p.value, "fdr")

# No significant indicative species 


### ______######
##### Tax Glom #####
# load("ps.merg.tree.RData")
# 
# ps_sp_glom <- ps.merg.tree
# 
# tax <- as.data.frame(tax_table(ps_sp_glom))
# 
# for (k in (1:7)){
#   for (i in 1:nrow(tax)){
#     if ((k=1) && is.na(tax[i,k])){
#       tax[i,k] <- paste0("K_NA_",i)
#     } else if ((k=2) && is.na(tax[i,k])){
#       tax[i,k] <- paste0("P_NA_",i)
#     } else if  ((k=3) && is.na(tax[i,k])){
#       tax[i,k] <- paste0("C_NA_",i)
#     } else if  ((k=4) && is.na(tax[i,k])){
#       tax[i,k] <- paste0("O_NA_",i)
#     } else if ((k=5) && is.na(tax[i,k])){
#       tax[i,k] <- paste0("F_NA_",i)
#     } else if  ((k=6) && is.na(tax[i,k])){
#       tax[i,k] <- paste0("G_NA_",i)
#     } else if  ((k=7) && is.na(tax[i,k])){
#       tax[i,k] <- paste0("S_NA_",i)}
#   }
# }
# 
# TAX <- tax_table(as.matrix(tax))
# tax_table(ps_sp_glom) <-TAX
# tax2 <- as.data.frame(tax_table(ps_sp_glom))
# identical(tax,tax2)
# 
# ps_ge_glom <-tax_glom(ps_sp_glom, taxrank="Genus",NArm=F) 
# save(ps_ge_glom,file="ps_ge_glom.RData")
# 
# ps_sp_glom <- tax_glom(ps_sp_glom, taxrank="Species",NArm=F)
# save(ps_sp_glom,file="ps_sp_glom.RData")

load("ps_sp_glom.RData")
ASV_ab_glom <- as.data.frame(otu_table(ps_sp_glom))
Tax_tab_glom <- as.data.frame(tax_table(ps_sp_glom))
Sdata <- as.data.frame(sample_data(ps_sp_glom))


ASV_fun_indval_glom <- multipatt(ASV_ab_glom,Sdata$Burning_intensity,
                            control=how(nperm=999), 
                            duleg =TRUE,
                            print.perm =TRUE)

summary(ASV_fun_indval_glom)
P_adj_glom <- ASV_fun_indval_glom$sign
P_adj_glom$p_val_adj_holm <- p.adjust(P_adj_glom$p.value, "holm")
P_adj_glom$p_val_adj_fdr <- p.adjust(P_adj_glom$p.value, "fdr")

signif <- subset(P_adj_glom,P_adj_glom$p.value<0.05)

signif_tax <- merge(signif,Tax_tab_glom,by=0)
row.names(signif_tax) <- signif_tax$Row.names

signif_tax_ab <- merge(signif_tax,t(ASV_ab_glom),by=0)

# 4 ASV ind of MedBurn
# 7 ASV ind of HighBurn


ASV_fun_indval_comb_glom <- multipatt(ASV_ab_glom,Sdata$Burning_intensity,
                                 control=how(nperm=999), 
                                 duleg =FALSE,
                                 restcomb = c(1,2,3,4,5,8,9,10,14), # Take only relevant combinations 
                                 print.perm =TRUE)

P_adj_comb_glom <- ASV_fun_indval_comb_glom$sign
P_adj_comb_glom$p_val_adj_holm <- p.adjust(P_adj_comb_glom$p.value, "holm")
P_adj_comb_glom$p_val_adj_fdr <- p.adjust(P_adj_comb_glom$p.value, "fdr")


signif_comb <- subset(P_adj_comb_glom,P_adj_comb_glom$p.value<0.05)
signif_tax_comb <- merge(signif_comb,Tax_tab_glom,by=0)
row.names(signif_tax_comb) <- signif_tax_comb$Row.names

signif_tax_ab_comb <- merge(signif_tax_comb,t(ASV_ab_glom),by=0)
row.names(signif_tax_ab_comb) <- signif_tax_ab_comb$Row.names
signif_tax_ab_comb <- signif_tax_ab_comb[,-c(1,2)]
# 10 ind fungi


summary(ASV_fun_indval_comb_glom,indvalcomp = T)

# 3 ind of High Burn :
# p val = 0.04 -> A = 0.99 / B = 1 
# p val = 0.037 -> A = 0.8455 / B = 1 
# p val = 0.0316 -> A = 0.9032 / B = 0.875

# 1 ind of no - low burn -> Only present in this group A= 1 / B + 0/75 
# Herpotrcichillaceae family

# 1 ind of low-med burn
# p val = 0.044 -> A = 0.862 / B = 0.875 


# 4 ind of med -high burn

# p val = 0.004 -> A = 0.99 / B = 1 
# p val = 0.037 -> A = 0.8455 / B = 1 
# p val = 0.016 -> A = 0.9032 / B = 0.875 
# p val = 0.015 -> A = 1 / B = 0.75 

# 1 ind of burn 
# p val = 0.022 -> A = 0.9427 / B = 1 



### ______######
### Filtering ASV > 0.1% ######
load("ps_sp_glom.RData")
minTotRelAbun = 1e-3
x = taxa_sums(ps_sp_glom)
keepTaxa = taxa_names(ps_sp_glom)[which((x / sum(x)) > minTotRelAbun)]
ps_sp_glom0.1pc = prune_taxa(keepTaxa, ps_sp_glom)
# Keeps 114 ASVs

ASV_0.1_ab <- as.data.frame(otu_table(ps_sp_glom0.1pc))
Tax_0.1_tab <- as.data.frame(tax_table(ps_sp_glom0.1pc))
Sdata_0.1 <- as.data.frame(sample_data(ps_sp_glom0.1pc))

# library(labdsv)
# ASV_bac_indval<-indval(ASV_ab,Sdata$Burning_intensity)
# summary(ASV_bac_indval)

ASV_fun_indval_0.1 <- multipatt(ASV_0.1_ab,Sdata_0.1$Burning_intensity,
                            control=how(nperm=999), 
                            duleg =TRUE,
                            print.perm =TRUE)

summary(ASV_fun_indval_0.1)
P_adj_0.1 <- ASV_fun_indval_0.1$sign
P_adj_0.1 $p_val_adj_holm <- p.adjust(P_adj_0.1 $p.value, "holm")
P_adj_0.1 $p_val_adj_fdr <- p.adjust(P_adj_0.1 $p.value, "fdr")


# No significant indicative specie for FDR 



ASV_fun_indval_0.1_comb <- multipatt(ASV_0.1_ab,Sdata_0.1$Burning_intensity,
                                      control=how(nperm=999), 
                                      duleg =FALSE,
                                      restcomb = c(1,2,3,4,5,8,9,10,14), # Take only combinations that make sense 
                                      print.perm =TRUE)

P_adj_0.1_comb <- ASV_fun_indval_0.1_comb$sign
P_adj_0.1_comb $p_val_adj_holm <- p.adjust(P_adj_0.1 $p.value, "holm")
P_adj_0.1_comb $p_val_adj_fdr <- p.adjust(P_adj_0.1 $p.value, "fdr")


### Filtering ASV > 0.5% ######
load("ps_sp_glom.RData")
minTotRelAbun = 5e-3
x = taxa_sums(ps_sp_glom)
keepTaxa = taxa_names(ps_sp_glom)[which((x / sum(x)) > minTotRelAbun)]
ps_sp_glom0.5pc = prune_taxa(keepTaxa, ps_sp_glom)
# Keeps 33 ASVs

ASV_0.5_ab <- as.data.frame(otu_table(ps_sp_glom0.5pc))
Tax_0.5_tab <- as.data.frame(tax_table(ps_sp_glom0.5pc))
Sdata_0.5 <- as.data.frame(sample_data(ps_sp_glom0.5pc))


ASV_fun_indval_0.5 <- multipatt(ASV_0.5_ab,Sdata_0.5$Burning_intensity,
                                control=how(nperm=999), 
                                duleg =TRUE,
                                print.perm =TRUE)

summary(ASV_fun_indval_0.5)
P_adj_0.5 <- ASV_fun_indval_0.5$sign
P_adj_0.5 $p_val_adj_holm <- p.adjust(P_adj_0.5 $p.value, "holm")
P_adj_0.5 $p_val_adj_fdr <- p.adjust(P_adj_0.5 $p.value, "fdr")

# No significant indicative specie for FDR 



ASV_fun_indval_0.5_comb <- multipatt(ASV_0.5_ab,Sdata_0.5$Burning_intensity,
                                     control=how(nperm=999), 
                                     duleg =FALSE,
                                     restcomb = c(1,2,3,4,5,8,9,10,14), # Take only combinations that make sense 
                                     print.perm =TRUE)

P_adj_0.5_comb <- ASV_fun_indval_0.5_comb$sign
P_adj_0.5_comb $p_val_adj_holm <- p.adjust(P_adj_0.5 $p.value, "holm")
P_adj_0.5_comb $p_val_adj_fdr <- p.adjust(P_adj_0.5 $p.value, "fdr")
# No significant indicative specie for FDR 

