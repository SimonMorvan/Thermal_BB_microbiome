#### GETTING READY ####
R.version #4.1.1 
rm(list = ls())
gc()

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#                            Phylogenetic tree script                                                                    
#                                                                            
#    This script aims to create a phylogenetic tree 
#    then to attach it to the phyloseq object 
#                                                                            
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#

# This script was run on a remote server. 


#### Making the tree #####

library(phyloseq)
library(dada2)
library(phangorn)
library(DECIPHER)

setwd(dir = "/storage/Simon_morvan/Reads/ch2/16S/Phylo_tree")
load("ps_bac_merg.RData")

ASV_seqs <- getSequences(otu_table(ps.bac.merg, taxa_are_rows=FALSE))
names(ASV_seqs) <- ASV_seqs


gT_seqs <- lapply(order(width(ASV_seqs), decreasing=TRUE), # Guidetree
                    function(x) {
                      attr(x, "height") <- 0
                      attr(x, "label") <- names(ASV_seqs)[x]
                      attr(x, "members") <- 1L
                      attr(x, "leaf") <- TRUE
                      x
                    })
attr(gT_seqs, "height") <- 0.5
attr(gT_seqs, "members") <- length(gT_seqs)
class(gT_seqs) <- "dendrogram"
head(gT_seqs)
save(gT_seqs, file = 'gT_seqs.RData')


# Alignment 
alignment <- AlignSeqs(DNAStringSet(ASV_seqs), 
                       anchor=NA,
                       guideTree=gT_seqs, 
                       processors = NULL, # Automatically detects processors
                       iterations=0,
                       refinements=0) 
Bac_align <- alignment
save(Bac_align, file = "Bac_align.RData")

# BrowseSeqs(ITS_align, highlight=0) #look at alignment
# alignment <- AdjustAlignment(ITS_align) #to make manual adjustments

phang.align <- phyDat(as(Bac_align, "matrix"), type="DNA")
#transforms alignment into phyDat format


dm <- dist.ml(phang.align) # Compute pairwise distance on the phyDat object 
Bac_dm <- dm
save(Bac_dm, file = "Bac_dm.RData")


treeNJ <- NJ(Bac_dm) # Neighbor-Joining tree
# Note, tip order != sequence order
save(treeNJ, file = 'Bac_treeNJ.RData')



fit <- pml(treeNJ, data = phang.align, 
           multicore = TRUE, mc.cores = NULL)
#calculate the likelihood for the tree given the data


fitGTR <- update(fit, k = 4, inv = 0.2)


fitGTR <- optim.pml(fitGTR, model = "GTR", 
                    optInv = TRUE, optGamma = TRUE,
                    rearrangement = "stochastic",
                    control = pml.control(trace = 0))
# optim.pml optimizes the different model parameters


Bac_NJ_GTR <- fitGTR
save(Bac_NJ_GTR, file = 'Bac_NJ_GTR.RData')

set.seed(11)
Bac_bs <- bootstrap.pml(Bac_NJ_GTR,  
                           bs = 100, #number of bs
                           optNni = TRUE, #optimize topology
                           trees = TRUE, #returns trees
                           multicore = TRUE, mc.cores = NULL) 

save(Bac_bs, file = 'Bac_bs.RData')
 

#### ON LOCAL COMPUTER #### 


# First download ITS_NJ_GTR to computer and load it
setwd(dir = "/Users/simonmorvan/Desktop/UdeM/Projet/Ch2/16S_Data_Analysis/03_RData/Phylo_tree/")
load("Bac_NJ_GTR.RData")


# Add the tree to ps object merged data
setwd(dir = "/Users/simonmorvan/Desktop/UdeM/Projet/Ch2/16S_Data_Analysis/03_RData/")
load("ps_bac_merg.RData")
ps.bacmerg.tree <- merge_phyloseq(ps.bac.merg,phy_tree(Bac_NJ_GTR$tree))
save(ps.bacmerg.tree, file = 'ps.bacmerg.tree.RData')

# Add the tree to ps object unmerged data
setwd(dir = "/Users/simonmorvan/Desktop/UdeM/Projet/Ch2/16S_Data_Analysis/03_RData/")
load("ps_bac_unmerged.RData")
ps.bacunmerged.tree <- merge_phyloseq(ps.bac.samples,phy_tree(Bac_NJ_GTR$tree))
save(ps.bacunmerged.tree, file = 'ps.bacunmerged.tree.RData')

###  Check if the tree make sense #####

load("ps.bacmerg.tree.RData")
ntaxa(ps.bacmerg.tree)
ps_rhizo <- subset_taxa(ps.bacmerg.tree, Order=="Rhizobiales")
tax <- as.data.frame(tax_table(ps_rhizo))

ps_Beijerinckiaceae <- subset_taxa(ps.bacmerg.tree, Family=="Beijerinckiaceae")
plot_tree(ps_Beijerinckiaceae, shape="Burning_intensity", color="Genus",label.tips="Species")

