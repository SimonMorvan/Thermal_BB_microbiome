#### GETTING READY ####
R.version #4.1.1 
rm(list = ls())
gc()

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#                           Metacoder script                                 #                                   
#                                                                            #                   
#                                                                            #
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#

setwd("/Users/simonmorvan/Desktop/UdeM/Projet/Ch2/16S_Data_Analysis/03_RData/")

##### Packages #####
library(phyloseq)
library(metacoder)

##### Data #####
load("ps_bac_merg.RData")


#### Metacoder object  ####
ps.bac.merg.relab <- transform_sample_counts(ps.bac.merg,function(x) x/sum(x))
# First convert abundance to relative abundance

meta_obj <- parse_phyloseq(ps.bac.merg.relab)
meta_obj$data$sample_data$Burning_intensity <- as.factor(meta_obj$data$sample_data$Burning_intensity)
meta_obj$data$sample_data$Plot <- as.factor(meta_obj$data$sample_data$Plot)



#### Differential heat tree #####


meta_obj$data$tax_abund <- calc_taxon_abund(meta_obj, "otu_table")
# Quantifies each taxon abundance -> sums the ASVs abundance with the same taxonomy


meta_obj$data$diff_table <- compare_groups(meta_obj, data = "tax_abund",
                                      cols = meta_obj$data$sample_data$Plot,# What columns of sample data to use
                                      groups = meta_obj$data$sample_data$Burning_intensity) # What category each sample is assigned to
# Compares the taxonomic difference based on ASV abundance

meta_obj$data$diff_table$wilcox_p_value_fdr <- p.adjust(meta_obj$data$diff_table$wilcox_p_value,
                                               method = "fdr")
# All the corrected p_values are not significant

print(meta_obj$data$diff_table$wilcox_p_value_fdr )
hist(meta_obj$data$diff_table$wilcox_p_value_fdr )

meta_obj$data$diff_table$log2_median_ratio[meta_obj$data$diff_table$wilcox_p_value_fdr > 0.1] <- 0  #start at 0.1
print(meta_obj$data$diff_table$log2_median_ratio)


heat_tree_matrix(meta_obj,
                 data = "diff_table",
                 node_size = n_obs, # n_obs is a function that calculates, in this case, the number of ASV per taxon
                 node_label = taxon_names,
                 node_color = log2_median_ratio, # A column from `obj$data$diff_table`
                 node_color_range = diverging_palette(), # The built-in palette for diverging data
                 node_color_trans = "linear", # The default is scaled by circle area
                 node_color_interval = c(-3, 3), # The range of `log2_median_ratio` to display
                 edge_color_interval = c(-3, 3), # The range of `log2_median_ratio` to display
                 node_size_axis_label = "Number of ASVs",
                 node_color_axis_label = "Log2 ratio median proportions",
                 layout = "davidson-harel", # The primary layout algorithm
                 initial_layout = "reingold-tilford", # The layout algorithm that initializes node locations
                 output_file = "differential_heat_tree.pdf") # Saves the plot as a pdf file



#### Complete heat tree #####
meta_obj <- parse_phyloseq(ps.bac.merg)
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
          title = "Bacterial dataset",
          output_file="16S_heat_tree.pdf")

### Check ###
Taxo <- as.data.frame(meta_obj$data$tax_data)
Tax_abu <- as.data.frame(meta_obj$data$tax_abund)


#### Relative ab. of specific taxa ####
bac_ab <- rowSums(otu_table(ps.bac.merg))
tot_bac_ab <- sum(bac_ab)

### Actinobacteriota ####
ps.act <- subset_taxa(ps.bac.merg,Phylum=="Actinobacteriota")
act_ab <- rowSums(otu_table(ps.act))
tot_act_ab <- sum(act_ab)
tot_act_ab/tot_bac_ab #  21.98

# Frankiales
ps.frankiales <- subset_taxa(ps.bac.merg,Order=="Frankiales")
frank_ab <- rowSums(otu_table(ps.frankiales))
tot_frank_ab <- sum(frank_ab)
tot_frank_ab/tot_bac_ab # -> 14.00% 

# Acidothermus
ps.acidothermus <- subset_taxa(ps.bac.merg,Genus=="Acidothermus")
actherm_ab <- rowSums(otu_table(ps.acidothermus))
tot_actherm_ab<- sum(actherm_ab)

bac_ab <- rowSums(otu_table(ps.bac.merg))
tot_bac_ab <- sum(bac_ab)

tot_actherm_ab/tot_bac_ab # -> 12.86% 

### Acidobacteria ####

ps.ac <- subset_taxa(ps.bac.merg,Phylum=="Acidobacteriota")
ac_ab <- rowSums(otu_table(ps.ac))
tot_ac_ab <- sum(ac_ab)
tot_ac_ab/tot_bac_ab #  21.98

# Acidobacteriales -> acz
# acz -> 12.56 % 

ps.acidobac <- subset_taxa(ps.bac.merg,Order=="Acidobacteriales")
acidobac_ab <- rowSums(otu_table(ps.acidobac))
tot_acidobac_ab <- sum(acidobac_ab)

tot_acidobac_ab/tot_bac_ab # -> 12.48% 

### Planctomycetota ####
ps.plancto <- subset_taxa(ps.bac.merg,Phylum=="Planctomycetota")
plancto_ab <- rowSums(otu_table(ps.plancto))
tot_plancto_ab <- sum(plancto_ab)
tot_plancto_ab/tot_bac_ab #  

# Isosphaerales 
# 12.7%
ps.iso <- subset_taxa(ps.bac.merg,Order=="Isosphaerales")
iso_ab <- rowSums(otu_table(ps.iso))
tot_iso_ab <- sum(iso_ab)
tot_iso_ab/tot_bac_ab # -> 12.41% 

ps.aqua <- subset_taxa(ps.bac.merg,Genus=="Aquisphaera")
aqua_ab <- rowSums(otu_table(ps.aqua))
tot_aqua_ab <- sum(aqua_ab)
tot_aqua_ab/tot_bac_ab # -> 10.1%

### Proteobacteria ####
ps.proteo <- subset_taxa(ps.bac.merg,Phylum=="Proteobacteria")
proteo_ab <- rowSums(otu_table(ps.proteo))
tot_proteo_ab <- sum(proteo_ab)
tot_proteo_ab/tot_bac_ab # -> 18.91% 

# Rhizobiales 
ps.rhizo <- subset_taxa(ps.bac.merg,Order=="Rhizobiales")
rhizo_ab <- rowSums(otu_table(ps.rhizo))
tot_rhizo_ab <- sum(rhizo_ab)
tot_rhizo_ab/tot_bac_ab # -> 12.41% 

# Roseiarcus 
ps.rose <- subset_taxa(ps.bac.merg,Genus=="Roseiarcus")
rose_ab <- rowSums(otu_table(ps.rose))
tot_rose_ab <- sum(rose_ab)
tot_rose_ab/tot_bac_ab # -> 10.1%

# Gemmatales
ps.gemma <- subset_taxa(ps.bac.merg,Order=="Gemmatales")
gemma_ab <- rowSums(otu_table(ps.gemma))
tot_gemma_ab <- sum(gemma_ab)
tot_gemma_ab/tot_bac_ab # -> 2.65% 

# Acidobacteriaceae (Subgroup 1)
ps.Subgroup1 <- subset_taxa(ps.bac.merg,Family=="Acidobacteriaceae (Subgroup 1)")
Subgroup1_ab <- rowSums(otu_table(ps.Subgroup1))
tot_Subgroup1_ab <- sum(Subgroup1_ab)
tot_Subgroup1_ab/tot_bac_ab # -> 6.2%

# Xanthobacteraceae 
ps.Xanthobacteraceae <- subset_taxa(ps.bac.merg,Family=="Xanthobacteraceae")
Xanthobacteraceae_ab <- rowSums(otu_table(ps.Xanthobacteraceae))
tot_Xanthobacteraceae_ab <- sum(Xanthobacteraceae_ab)
tot_Xanthobacteraceae_ab/tot_bac_ab # 


