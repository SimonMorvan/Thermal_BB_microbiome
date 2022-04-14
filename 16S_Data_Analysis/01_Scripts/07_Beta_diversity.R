# Simon Morvan
# Sept 2O21
R.version # 4.1.1.

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#                                                                          # 
#                               Beta diversity                             # 
#                                                                          #
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
Dir <- "/Users/simonmorvan/Desktop/UdeM/Projet/Ch2/16S_Data_Analysis/03_RData/"
setwd(Dir)

#### Packages ####
library(phyloseq)
library(ggplot2)
library(ggpubr)
library(vegan)
library(microbiome)
library(ggrepel)

#### Data ####
load("ps.bacmerg.tree.RData")

###########______________##########
####CLR#####
# Based on 
# https://www.nicholas-ollberding.com/post/introduction-to-the-statistical-analysis-of-microbiome-data-in-r/


# Centered log-ratio transformation
ps_clr <- microbiome::transform(ps.bacmerg.tree, "clr")        

#### PCA ####
ord_clr <- phyloseq::ordinate(ps_clr, "RDA", distance = "euclidean")

phyloseq::plot_scree(ord_clr) + 
  geom_bar(stat="identity", fill = "blue") +
  labs(x = "\nAxis", y = "Proportion of Variance\n")

head(ord_clr$CA$eig) 
sapply(ord_clr$CA$eig[1:5], function(x) x / sum(ord_clr$CA$eig))     

#Scale axes and plot ordination
clr1 <- ord_clr$CA$eig[1] / sum(ord_clr$CA$eig)
clr2 <- ord_clr$CA$eig[2] / sum(ord_clr$CA$eig)

AitBetaDivBac<- phyloseq::plot_ordination(ps_clr, ord_clr, type="samples", color="Burning_intensity") + 
                theme_bw()+
                ggtitle(label = "Aitchison distance - 16S")+
                geom_point(colour = "black", size = 4.5) +
                geom_point(size = 4)+
                geom_text_repel(aes(label=Plot), color="black")+
                theme(axis.text=element_text(size=12))+ # size of lables, axis titles
                theme(axis.title.x = element_text(size=15), axis.title.y = element_text(size=15))+
                theme(legend.text=element_text(size=15),
                      legend.title=element_text(size=15),
                      legend.background = element_rect(size=0.5, linetype="solid", 
                                                       colour ="black"))+
                scale_color_manual(values = c("#9cccd8","#ffb53f","#f26419","#931125"))


#### PERMANOVA #### 

#Generate distance matrix
clr_dist_matrix <- phyloseq::distance(ps_clr, method = "euclidean") 

metadata <- as(sample_data(ps_clr), "data.frame") # Export data

perm <- how(nperm = 999)
setBlocks(perm)<-with(metadata,Block)
Permanova.clr<-adonis2(clr_dist_matrix~Burning_intensity,
                       data = metadata,
                       permutations = perm)


#Dispersion test and plot
betadisp <- betadisper(clr_dist_matrix, metadata$Burning_intensity)
permutest(betadisp, permutations = perm) 

plot(betadisp, main = "Ordination Centroids and Dispersion Labeled: Aitchison Distance", sub = "")
boxplot(betadisp, main = "", xlab = "")


###########______________##########
#### Hellinger ####
sum(otu_table(ps.bacmerg.tree)==0) # 8952  zeros
sum(otu_table(ps.bacmerg.tree)==0)/(nrow(otu_table(ps.bacmerg.tree))*ncol(otu_table(ps.bacmerg.tree)))*100 
# 59 % of the dataframe is composed of zeros

ps.bac.hell <- transform_sample_counts(ps.bacmerg.tree, function(x) sqrt(x/sum(x)))
# Hellinger transformation

# Check
# hell.2<- as.data.frame(otu_table(ps.bac.hell))
# Te <- (as.data.frame(otu_table(ps.bac.merg)))
# hell <- as.data.frame(decostand(Te, method="hellinger"))# ASV as columns! 

#### PCA ####
ord.PCA.hellinger <- ordinate(ps.bac.hell, method="RDA",distance = "euclidean")

phyloseq::plot_scree(ord.PCA.hellinger) + 
  geom_bar(stat="identity", fill = "blue") +
  labs(x = "\nAxis", y = "Proportion of Variance\n")

PCA_Hell_bac <- plot_ordination(ps.bac.hell, ord.PCA.hellinger, type="samples",color = "Burning_intensity")+ 
  theme_bw()+
  ggtitle(label = "Hellinger distance - 16S")+
  geom_point(colour = "black", size = 4.5) +
  geom_point(size = 4)+
  geom_text_repel(aes(label=Plot), color="black")+
  theme(axis.text=element_text(size=12))+ # size of lables, axis titles
  theme(axis.title.x = element_text(size=15), axis.title.y = element_text(size=15))+
  theme(legend.text=element_text(size=15),
        legend.title=element_text(size=15),
        legend.background = element_rect(size=0.5, linetype="solid", 
                                         colour ="black"))+
  scale_color_manual(values = c("#9cccd8","#ffb53f","#f26419","#931125"))
PCA_Hell_bac


#### PERMANOVA ####
metadata <- as(sample_data(ps.bac.hell), "data.frame") # Export data

perm <- how(nperm = 999)
setBlocks(perm)<-with(metadata,Block)
Permanova.hell <-adonis2(phyloseq::distance(physeq = ps.bac.hell, method="euclidean")~Burning_intensity,
                         data = metadata,
                         permutations = perm)



betadisp <- betadisper(phyloseq::distance(physeq = ps.bac.hell, method="euclidean"), metadata$Burning_intensity)
permutest(betadisp, permutations = perm) 

# The betadisp tests homogeneity of dispersion among groups (regions in your case), 
# which is a condition (assumption) for adonis. 
# You want the permutest to be not significant 
# because the null hypothesis is that there is no difference in terms of disperison

# You may have the centroids of two groups in NMS at a very similar position in the ordination space, 
# but if their dispersions are quite different, adonis will give you a significant p-value, 
# thus, the result is heavily influenced not by the difference in composition 
# between groups but by differences in composition within groups (heterogeneous dispersion and thus a measure of betadiversity). 
# In short, your results are fine, you are meeting the 'one assumption' for adonis (homogeneous dispersion)
# and thus you are certain that results from adonis are 'real' and not an artifact of heterogeneous dispersions. 

###########______________##########
#### Unweighted UNIFRAC ####

#### PCoA ####
ord.PCoA.Unifrac <- ordinate(ps.bacmerg.tree, method="PCoA", distance = "uunifrac")

PCoA_Unifrac_BAC <- plot_ordination(ps.bacmerg.tree, ord.PCoA.Unifrac, type="samples",color = "Burning_intensity")+ 
  theme_bw()+
  ggtitle(label = "Unweighted Unifrac distance - 16S")+
  geom_point(colour = "black", size = 4.5) +
  geom_point(size = 4)+
  geom_text_repel(aes(label=Plot), color="black")+
  theme(axis.text=element_text(size=12))+ # size of lables, axis titles
  theme(axis.title.x = element_text(size=15), axis.title.y = element_text(size=15))+
  theme(legend.text=element_text(size=15),
        legend.title=element_text(size=15),
        legend.background = element_rect(size=0.5, linetype="solid", 
                                         colour ="black"))+
  scale_color_manual(values = c("#9cccd8","#ffb53f","#f26419","#931125"))
PCoA_Unifrac_BAC


#### PERMANOVA ####
metadata <- as(sample_data(ps.bacmerg.tree), "data.frame") # Export data

perm <- how(nperm = 999)
setBlocks(perm)<-with(metadata,Block)
Permanova.Uni <-adonis2(phyloseq::distance(physeq = ps.bacmerg.tree, method="uunifrac")~Burning_intensity,
                         data = metadata,
                         permutations = perm)



betadisp <- betadisper(phyloseq::distance(physeq = ps.bacmerg.tree, method="uunifrac"), metadata$Burning_intensity)
permutest(betadisp, permutations = perm) 


###########______________##########

#### Weighted UNIFRAC ####

#### PCoA ####
ord.PCoA.WUnifrac <- ordinate(ps.bacmerg.tree, method="PCoA", distance = "wunifrac")

PCoA_WUnifrac_BAC <- plot_ordination(ps.bacmerg.tree, ord.PCoA.WUnifrac, type="samples",color = "Burning_intensity")+ 
  theme_bw()+
  ggtitle(label = "Weighted Unifrac distance - 16S")+
  geom_point(colour = "black", size = 4.5) +
  geom_point(size = 4)+
  geom_text_repel(aes(label=Plot), color="black")+
  theme(axis.text=element_text(size=12))+ # size of lables, axis titles
  theme(axis.title.x = element_text(size=15), axis.title.y = element_text(size=15))+
  theme(legend.text=element_text(size=15),
        legend.title=element_text(size=15),
        legend.background = element_rect(size=0.5, linetype="solid", 
                                         colour ="black"))+
  scale_color_manual(values = c("#9cccd8","#ffb53f","#f26419","#931125"))
PCoA_WUnifrac_BAC



#### PERMANOVA ####
metadata <- as(sample_data(ps.bacmerg.tree), "data.frame") # Export data

perm <- how(nperm = 999)
setBlocks(perm)<-with(metadata,Block)
Permanova.WUni <-adonis2(phyloseq::distance(physeq = ps.bacmerg.tree, method="wunifrac")~Burning_intensity,
                        data = metadata,
                        permutations = perm)



betadisp <- betadisper(phyloseq::distance(physeq = ps.merg.tree, method="wunifrac"), metadata$Burning_intensity)
permutest(betadisp, permutations = perm) 


###########______________##########

#### Bray Curtis - log #####

ps.bac.log <- transform_sample_counts(ps.bacmerg.tree, function(x) log(x+1))

###NMDS####
ord.nmds.bray <- ordinate(ps.bac.log, method="NMDS",  distance = "bray", k = 2, try = 100) 

NMDS_Bray_log_BAC <- plot_ordination(ps.bac.log, ord.nmds.bray, type="samples",color = "Burning_intensity")+ 
  annotate("text", x=0.1, y=0.2 , label= "Stress = 0.17", size=5)+
  theme_bw()+
  ggtitle(label = "NMDS - 16S dataset _ Bray-curtis on log transformed data")+
  geom_point(colour = "black", size = 4.5) +
  geom_point(size = 4)+
  #geom_line()+
  geom_text_repel(aes(label=Plot), color="black")+
  theme(axis.text=element_text(size=12))+ # size of lables, axis titles
  theme(axis.title.x = element_text(size=15), axis.title.y = element_text(size=15))+
  theme(legend.text=element_text(size=15),
        legend.title=element_text(size=15),
        legend.background = element_rect(size=0.5, linetype="solid", 
                                         colour ="black"))+
  scale_color_manual(values = c("#9cccd8","#ffb53f","#f26419","#931125"))
NMDS_Bray_log_BAC

#### PERMANOVA ####
metadata <- as(sample_data(ps.bac.log), "data.frame") # Export data


perm <- how(nperm = 999)
setBlocks(perm)<-with(metadata,Block)
Permanova.log <-adonis2(phyloseq::distance(physeq = ps.bac.log, method="bray")~Burning_intensity,
                         data = metadata,
                         permutations = perm)


betadisp <- betadisper(phyloseq::distance(physeq = ps.bac.log, method="bray"), metadata$Burning_intensity)
permutest(betadisp, permutations = perm) 
