#### GETTING READY ####
R.version #4.1.1 
rm(list = ls())
gc()

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#                       Indicative species  script                           #                                   
#                                                                            #                   
#                                                                            #
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#


setwd("/Users/simonmorvan/Desktop/UdeM/Projet/Ch2/16S_Data_Analysis/03_RData/")

##### Packages #####
library(phyloseq)
library(vegan)
library(indicspecies)

##### Data #####
load("ps_bac_merg.RData")

ASV_ab <- as.data.frame(otu_table(ps.bac.merg))
Tax_tab <- as.data.frame(tax_table(ps.bac.merg))
Sdata <- as.data.frame(sample_data(ps.bac.merg))


ASV_bac_indval<- multipatt(ASV_ab,Sdata$Burning_intensity,
                               control=how(nperm=9999), 
                               duleg =TRUE,
                               print.perm =TRUE)

summary(ASV_bac_indval_comb)
P_adj <- ASV_bac_indval$sign
P_adj$p_val_adj_holm <- p.adjust(P_adj$p.value, "holm")
P_adj$p_val_adj_fdr <- p.adjust(P_adj$p.value, "fdr")

# No significant indicative species 

##### MULTIPATT comb sites ####

ASV_bac_indval_comb <- multipatt(ASV_ab,Sdata$Burning_intensity,
                            control=how(nperm=999), 
                            duleg =FALSE,
                            print.perm =TRUE)

P_adj_comb<- ASV_bac_indval_comb$sign
P_adj_comb$p_val_adj_holm <- p.adjust(P_adj_comb$p.value, "holm")
P_adj_comb$p_val_adj_fdr <- p.adjust(P_adj_comb$p.value, "fdr")

# No significant indicative species 



##### Tax Glom #####
#  load("ps_bac_merg.RData")
# # 
#  ps.sp.glom <- ps.bac.merg
# # 
# tax <- as.data.frame(tax_table(ps.sp.glom))
# tax[ tax == "NA" ] <- NA 
# 
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
# tax_table(ps.sp.glom) <-TAX
# tax2 <- as.data.frame(tax_table(ps.sp.glom))
# identical(tax,tax2)
# 
# ps.ge.glom <-tax_glom(ps.sp.glom, taxrank="Genus",NArm=F) 
# save(ps.ge.glom,file="ps.ge.glom.RData")
# 
# ps.sp.glom <- tax_glom(ps.sp.glom, taxrank="Species",NArm=F)
# save(ps.sp.glom,file="ps.sp.glom.RData")

load("ps.sp.glom.RData")
ASV_ab_glom <- as.data.frame(otu_table(ps.sp.glom))
Tax_tab_glom <- as.data.frame(tax_table(ps.sp.glom))
Sdata <- as.data.frame(sample_data(ps.sp.glom))


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

# NO ind spe

ASV_fun_indval_comb_glom <- multipatt(ASV_ab_glom,Sdata$Burning_intensity,
                                      control=how(nperm=999), 
                                      duleg =FALSE,
                                      restcomb = c(1,2,3,4,5,8,9,10,14), # Take only combinations that make sense 
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

# NO ind spe

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
load("ps.sp.glom.RData")
minTotRelAbun = 1e-3
x = taxa_sums(ps.sp.glom)
keepTaxa = taxa_names(ps.sp.glom)[which((x / sum(x)) > minTotRelAbun)]
ps.sp.glom0.1pc = prune_taxa(keepTaxa, ps.sp.glom)
# Keeps 193 ASVs

ASV_0.1_ab <- as.data.frame(otu_table(ps.sp.glom0.1pc))
Tax_0.1_tab <- as.data.frame(tax_table(ps.sp.glom0.1pc))
Sdata_0.1 <- as.data.frame(sample_data(ps.sp.glom0.1pc))


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
# No significant indicative specie for FDR 


### Filtering ASV > 0.5% ######
load("ps.sp.glom.RData")
minTotRelAbun = 5e-3
x = taxa_sums(ps.sp.glom)
keepTaxa = taxa_names(ps.sp.glom)[which((x / sum(x)) > minTotRelAbun)]
ps_sp_glom0.5pc = prune_taxa(keepTaxa, ps.sp.glom)
# Keeps 25 ASVs

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

