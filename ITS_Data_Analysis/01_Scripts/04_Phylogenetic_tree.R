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

setwd(dir = "/storage/Simon_morvan/Reads/ch2/ITS/Phylo_tree")
load("ps.merg.RData")
ASV_seqs <- getSequences(otu_table(ps.fun.merg, taxa_are_rows=FALSE))
names(ASV_seqs) <- ASV_seqs

gT_seqs <- lapply(order(width(ASV_seqs), decreasing=TRUE), #Guidetree
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



alignment <- AlignSeqs(DNAStringSet(seqs), 
                       anchor=NA,
                       guideTree=gT_seqs,
                       processors = NULL, # Automaticly detects processors
                       iterations=0,
                       refinements=0) 
ITS_align <- alignment
save(ITS_align, file = "ITS_align.RData")

# BrowseSeqs(ITS_align, highlight=0) #look at alignment
# alignment <- AdjustAlignment(ITS_align) #to make manual adjustments

phang.align <- phyDat(as(ITS_align, "matrix"), type="DNA")
#transforms alignment into phyDat format


dm <- dist.ml(phang.align) # Compute pairwise distance on the phyDat object 
ITS_dm <- dm
save(ITS_dm, file = "ITS_dm.RData")


treeNJ <- NJ(dm) # Neighbor-Joining tree
# Note, tip order != sequence order
save(treeNJ, file = 'ITS_treeNJ.RData')


fit <- pml(treeNJ, data = phang.align, 
           multicore = TRUE, mc.cores = NULL)
#calculate the likelihood for the tree given the data


fitGTR <- update(fit, k = 4, inv = 0.2)



fitGTR <- optim.pml(fitGTR, model = "GTR", 
                    optInv = TRUE, optGamma = TRUE,
                    rearrangement = "stochastic",
                    control = pml.control(trace = 0))
# optim.pml optimizes the different model parameters

ITS_NJ_GTR <- fitGTR
save(ITS_NJ_GTR, file = 'ITS_NJ_GTR.RData')

set.seed(11)
ITS_bs <- bootstrap.pml(ITS_NJ_GTR,  
                           bs = 100, #number of bs
                           optNni = TRUE, #optimize topology
                           trees = TRUE, #returns trees
                           multicore = TRUE, mc.cores = NULL) 
save(ITS_bs, file = 'ITS_bs.RData')


### ON LOCAL COMPUTER ####
# First download ITS_NJ_GTR to computer and load it
setwd(dir = "/Users/simonmorvan/Desktop/UdeM/Projet/Ch2/ITS_Data_Analysis/03_RData/Phylo_tree/")
load("ITS_NJ_GTR.RData")


# Add the tree to  ps object merged data
setwd(dir = "/Users/simonmorvan/Desktop/UdeM/Projet/Ch2/ITS_Data_Analysis/03_RData/")
load("ps.merg.RData")
ps.merg.tree <- merge_phyloseq(ps.fun.merg,phy_tree(ITS_NJ_GTR$tree))
save(ps.merg.tree, file = 'ps.merg.tree.RData')

# Add the tree to ps object unmerged data
setwd(dir = "/Users/simonmorvan/Desktop/UdeM/Projet/Ch2/ITS_Data_Analysis/03_RData/")
load("ps_fun_unmerged.RData")
ps.unmerged.tree <- merge_phyloseq(ps.fun.samples,phy_tree(ITS_NJ_GTR$tree))
save(ps.unmerged.tree, file = 'ps.unmerged.tree.RData')

###  Check if the tree make sense #####
load("ps.merg.tree.RData")
ntaxa(ps.merg.tree)
ps_helo <- subset_taxa(ps.merg.tree, Order=="Helotiales")
plot_tree(ps_helo, shape="Burning_intensity", color="Genus",label.tips="Species")
tax <- as.data.frame(tax_table(ps_helo))

ps_Hyaloscyphaceae <- subset_taxa(ps_helo, Family=="Hyaloscyphaceae")
plot_tree(ps_Hyaloscyphaceae, shape="Burning_intensity", color="Genus",label.tips="Species")


ps_Helotiaceae <- subset_taxa(ps_helo, Family=="Helotiaceae")
plot_tree(ps_Helotiaceae, shape="Burning_intensity", color="Genus",label.tips="Species")

