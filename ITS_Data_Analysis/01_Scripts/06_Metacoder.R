#### GETTING READY ####
R.version #4.1.1 
rm(list = ls())
gc()

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#                           Metacoder script                                 #                                   
#                                                                            #                   
#                                                                            #
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#

setwd("/Users/simonmorvan/Desktop/UdeM/Projet/Ch2/ITS_Data_Analysis/03_RData/")

##### Packages #####
library(phyloseq)
library(metacoder)

##### Data #####
load("ps.merg.tree.RData")


#### Metacoder object  ####
ps.fun.merg.relab <- transform_sample_counts(ps.merg.tree,function(x) x/sum(x))
# First convert abundance to relative abundance

meta_obj <- parse_phyloseq(ps.fun.merg.relab)
meta_obj$data$sample_data$Burning_intensity <- as.factor(meta_obj$data$sample_data$Burning_intensity)
meta_obj$data$sample_data$Plot <- as.factor(meta_obj$data$sample_data$Plot)




#### Diff heat tree #####


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





#### Full heat tree #####


meta_obj <- parse_phyloseq(ps.merg.tree)
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
          title = "Fungal dataset",
          output_file="ITS_heat_tree.pdf")

?heat_tree

### Check ###
Taxo <- as.data.frame(meta_obj$data$tax_data)
Tax_abu <- as.data.frame(meta_obj$data$tax_abund)
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$

# helotiales --> abl
#  relab -> ~46.0 %


ps.helotiales <- subset_taxa(ps.merg.tree,Order=="Helotiales")
helo_ab <- rowSums(otu_table(ps.helotiales))
tot_helo_ab <- sum(helo_ab)

fun_ab <- rowSums(otu_table(ps.merg.tree))
tot_fun_ab <- sum(fun_ab)

tot_helo_ab/tot_fun_ab # -> 46.8% 

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$

# Pezoloma ericae -> arl
# arl -> 12.5% 

ps.peri <- subset_taxa(ps.merg.tree,Species=="ericae")
peri_ab <- rowSums(otu_table(ps.peri))
tot_peri_ab <- sum(peri_ab)

tot_peri_ab/tot_fun_ab # -> 12.2% 



#$$$$$$$$$$$$$$$$$$$$$$$$$$$$

# O. maius ->  ark
# ark -> 9% 
ps.omai <- subset_taxa(ps.merg.tree,Species=="maius")
omai_ab <- rowSums(otu_table(ps.omai))
tot_omai_ab <- sum(omai_ab)

tot_omai_ab/tot_fun_ab # -> 9% 


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$

# Herpotrichiellaceae ->  aee
# ark -> 28.4 % 

ps.herpo <- subset_taxa(ps.merg.tree,Family=="Herpotrichiellaceae")
herpo_ab <- rowSums(otu_table(ps.herpo))
tot_herpo_ab <- sum(herpo_ab)

tot_herpo_ab/tot_fun_ab # -> 9% 


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$

# Agaricales ->  abm
# abm -> 6 % 

ps.agar <- subset_taxa(ps.merg.tree,Order=="Agaricales")
agar_ab <- rowSums(otu_table(ps.agar))
tot_agar_ab <- sum(agar_ab)

tot_agar_ab/tot_fun_ab # -> 5.6% 


# C. sphagnicola ->  arm
# abm -> 4 % 

ps.Csphagn <- subset_taxa(ps.agar,Species=="sphagnicola")
Csphagn_ab <- rowSums(otu_table(ps.Csphagn))
tot_Csphagn_ab <- sum(Csphagn_ab)

tot_Csphagn_ab/tot_fun_ab # -> 3.8% 
# Seems to work

ps.fortinii <- subset_taxa(ps.merg.tree,Species=="fortinii")
fortinii_ab <- rowSums(otu_table(ps.fortinii))
tot_fortiniin_ab <- sum(fortinii_ab)

tot_fortiniin_ab/tot_fun_ab # -> 3.8% 


ps.Chaetothyriales  <- subset_taxa(ps.merg.tree,Order=="Chaetothyriales")
Chaetothyriales_ab <- rowSums(otu_table(ps.Chaetothyriales))
tot_Chaetothyriales_ab  <- sum(Chaetothyriales_ab)

tot_Chaetothyriales_ab/tot_fun_ab # -> 28.6%% 

ps.Herpotrichiellaceae  <- subset_taxa(ps.merg.tree,Family=="Herpotrichiellaceae")
Herpotrichiellaceae_ab <- rowSums(otu_table(ps.Herpotrichiellaceae))
tot_Herpotrichiellaceae_ab  <- sum(Herpotrichiellaceae_ab)

tot_Herpotrichiellaceae_ab/tot_fun_ab # -> 27.3%% 
