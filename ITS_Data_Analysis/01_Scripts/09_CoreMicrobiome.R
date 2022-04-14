R.version #4.1.1 
rm(list = ls())
gc()

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#                            Core Microbiome                                 #                                   
#                                                                            #
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
Dir <- "/Users/simonmorvan/Desktop/UdeM/Projet/Ch2/ITS_Data_Analysis/03_RData/"
setwd(Dir)

##### Packages #####
library(MicEco)

load("ps.merg.tree.RData")
?ps_venn


ITS_Venn_ASV <- ps_venn(ps.merg.tree,
                        group = "Burning_intensity",
                        labels = list(col = "gray12", cex=1.2),
                        weight=F,
                        quantities = list(type=c("counts"),col = "gray12"),
                        col = "black",
                        fill= c("#9cccd8","#ffb53f","#f26419","#931125"),
                        adjust_labels=T)
ITS_Venn_ASV

ITS_Venn_relAb <- ps_venn(ps.merg.tree,
                        group = "Burning_intensity",
                        labels = list(col = "gray12", cex=1.5),
                        weight=T,
                        quantities = list(type=c("percent"),col = "gray12"),
                        col = "black",
                        fill= c("#9cccd8","#ffb53f","#f26419","#931125"),
                        adjust_labels=T)
ITS_Venn_relAb


ITS_Venn <- ps_venn(ps.merg.tree,
                        group = "Burning_intensity",
                        labels = list(col = "gray12", cex=1.1),
                        weight=F,
                        quantities = list(labels = c("144\n(2%)",
                                                     "69\n(1%)",
                                                     "101\n(2%)",
                                                     "86\n(3%)",
                                                     "18\n(1%)",
                                                     "14\n(2%)",
                                                     "5\n(<1%)",
                                                     "10\n(1%)",
                                                     "5\n(<1%)",
                                                     "10\n(2%)",
                                                     "47\n(3%)",
                                                     "10\n(<1%)",
                                                     "9\n(2%)",
                                                     "18\n(1%)",
                                                     "222\n(81%)")),
                        col = "black",
                        fill= c("#9cccd8","#ffb53f","#f26419","#931125"),
                        adjust_labels=T)
ITS_Venn


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


#### Full heat tree #####
meta_obj <- parse_phyloseq(ps.core)
meta_obj$data$sample_data$Burning_intensity <- as.factor(meta_obj$data$sample_data$Burning_intensity)
meta_obj$data$sample_data$Plot <- as.factor(meta_obj$data$sample_data$Plot)
meta_obj$data$sample_data$Type <- "All"
meta_obj$data$sample_data$Type <- as.factor(meta_obj$data$sample_data$Type)

meta_obj$data$relab <- calc_obs_props(meta_obj,"otu_table")
meta_obj$data$tax_abund <- calc_taxon_abund(meta_obj, "relab",
                                            cols = meta_obj$data$sample_data$sample_id,
                                            groups = meta_obj$data$sample_data$Type)
# Quantifies each taxon abundance -> sums the ASVs abundance with the same taxonomy across the whole dataset

meta_obj$data$tax_abund$Relab <- meta_obj$data$tax_abund$All/16
# Need to divide by the number of samples 

heat_tree(meta_obj,
          node_label = taxon_names,
          node_color= n_obs,
          node_color_axis_label="Number of ASVs",
          node_color_range = c("#FDE725FF", "#21908CFF","#440154FF"),
          node_size = (meta_obj$data$tax_abund$Relab),
          node_size_axis_label = "Relative abundance",
          title = "Core microbiome fungal dataset",
          output_file="ITSCore_heat_tree.pdf")

?heat_tree




# test <- as.data.frame(otu_table(ps.merg.tree))
# test2 <- as.data.frame(colSums(test))
# tot <- sum(test2)
# test3 <- test2/tot
# 
# tax <- as.data.frame(tax_table(ps.merg.tree))
# 
# test4 <- cbind(test3,tax)
# library(ggplot2)
# 
# 
# 
# sequences <- Biostrings::DNAStringSet(taxa_names(ps.merg.tree))
# names(sequences) <- taxa_names(ps.merg.tree)
# ps.merge.tree <- merge_phyloseq(ps.merg.tree, sequences)
# 
# taxa_names(ps.merge.tree) <- paste0("ASV", seq(ntaxa(ps.merge.tree)))
# refseq(ps.merge.tree)
# 

