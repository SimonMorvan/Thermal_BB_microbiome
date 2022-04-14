# Simon Morvan
# Sept 2O21
R.version # 4.1.1.
Dir <- "/Users/simonmorvan/Desktop/UdeM/Projet/Ch2/ITS_Data_Analysis/03_RData/"
setwd(Dir)

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#


#       Script for visualizing pseudoreplicate resemblance via ordinations   #                     


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#

#### Packages ####
library(phyloseq)
library(ggplot2)
library(ggrepel)
library(ggpubr)


#### ßdiv - Pseudorep  ####
load(file ="ps.unmerged.tree.RData")


#### I Hellinger distance ####
ps.fun.hell <- transform_sample_counts(ps.unmerged.tree, function(x) sqrt(x/sum(x)))

# PCA #
ord.PCA.hellinger <- ordinate(ps.fun.hell, method="RDA", distance = "euclidean")


PCA_Hell <- plot_ordination(ps.fun.hell, ord.PCA.hellinger, type="samples",color="Plot",shape="Burning_intensity")+ 
  theme_bw()+
  ggtitle(label = "PCA - Hellinger distance - ITS ")+
  geom_point(colour = "black", size = 4.5) +
  geom_point(size = 4)+
  geom_text_repel(aes(label=Sample_ID), color="black")+
  theme(axis.text=element_text(size=12))+ # size of lables, axis titles
  theme(axis.title.x = element_text(size=15), axis.title.y = element_text(size=15))+
  theme(legend.text=element_text(size=15),
        legend.title=element_text(size=15),
        legend.background = element_rect(size=0.5, linetype="solid", 
                                         colour ="black"))+
  geom_line(aes(group = Plot))
PCA_Hell


PCA_Hell_B <- plot_ordination(ps.fun.hell, ord.PCA.hellinger, type="samples", color = "Burning_intensity")+ 
  theme_bw()+
  ggtitle(label = "PCA - Hellinger distance - ITS ")+
  geom_point(colour = "black", size = 4.5) +
  geom_point(size = 4)+
  geom_text_repel(aes(label=Sample_ID), color="black")+
  theme(axis.text=element_text(size=12))+ # size of lables, axis titles
  theme(axis.title.x = element_text(size=15), axis.title.y = element_text(size=15))+
  theme(legend.text=element_text(size=15),
        legend.title=element_text(size=15),
        legend.background = element_rect(size=0.5, linetype="solid", 
                                         colour ="black"))+
  scale_color_manual(values = c("#9cccd8","#ffb53f","#f26419","#931125"))+
  geom_line(aes(group = Plot))
PCA_Hell_B


#### I CLR distance ####

library(microbiome)


# Centered log-ratio transformation
ps_clr <- microbiome::transform(ps.unmerged.tree, "clr")        


#PCA #
ord_clr <- phyloseq::ordinate(ps_clr, "RDA", distance = "Plot")


phyloseq::plot_scree(ord_clr) + 
  geom_bar(stat="identity", fill = "blue") +
  labs(x = "\nAxis", y = "Proportion of Variance\n")

head(ord_clr$CA$eig) 
sapply(ord_clr$CA$eig[1:5], function(x) x / sum(ord_clr$CA$eig))     

#Scale axes and plot ordination
clr1 <- ord_clr$CA$eig[1] / sum(ord_clr$CA$eig)
clr2 <- ord_clr$CA$eig[2] / sum(ord_clr$CA$eig)

# PCA #

PCA_Ait <- plot_ordination(ps_clr, ord_clr, type="samples",color="Plot",shape="Burning_intensity") + 
  theme_bw()+
  ggtitle(label = "PCA - Aitchison distance - ITS")+
  geom_point(colour = "black", size = 4.5) +
  geom_point(size = 4)+
  geom_text_repel(aes(label=Sample_ID), color="black")+
  theme(axis.text=element_text(size=12))+ # size of lables, axis titles
  theme(axis.title.x = element_text(size=15), axis.title.y = element_text(size=15))+
  theme(legend.text=element_text(size=15),
        legend.title=element_text(size=15),
        legend.background = element_rect(size=0.5, linetype="solid", 
                                         colour ="black"))+
  geom_line(aes(group = Plot))
PCA_Ait


PCA_Ait_B <- plot_ordination(ps_clr, ord_clr, type="samples",color="Burning_intensity") + 
  theme_bw()+
  ggtitle(label = "PCA - Aitchison distance - ITS")+
  geom_point(colour = "black", size = 4.5) +
  geom_point(size = 4)+
  geom_text_repel(aes(label=Sample_ID), color="black")+
  theme(axis.text=element_text(size=12))+ # size of lables, axis titles
  theme(axis.title.x = element_text(size=15), axis.title.y = element_text(size=15))+
  theme(legend.text=element_text(size=15),
        legend.title=element_text(size=15),
        legend.background = element_rect(size=0.5, linetype="solid", 
                                         colour ="black"))+
  scale_color_manual(values = c("#9cccd8","#ffb53f","#f26419","#931125"))+
  geom_line(aes(group = Plot))
PCA_Ait_B


#### I Unifrac distance ####

ord.PCoA.Unifrac <- ordinate(ps.unmerged.tree, method="PCoA", distance = "uunifrac")

PCoA_Unifrac <- plot_ordination(ps.unmerged.tree, ord.PCoA.Unifrac, type="samples", color = "Plot",shape="Burning_intensity")+ 
  theme_bw()+
  ggtitle(label = "PCoA - Unweigthed Unifrac distance - ITS ")+
  geom_point(colour = "black", size = 4.5) +
  geom_point(size = 4)+
  geom_text_repel(aes(label=Sample_ID), color="black")+
  theme(axis.text=element_text(size=12))+ # size of lables, axis titles
  theme(axis.title.x = element_text(size=15), axis.title.y = element_text(size=15))+
  theme(legend.text=element_text(size=15),
        legend.title=element_text(size=15),
        legend.background = element_rect(size=0.5, linetype="solid", 
                                         colour ="black"))+
  geom_line(aes(group = Plot))
PCoA_Unifrac

PCoA_Unifrac_B <- plot_ordination(ps.unmerged.tree, ord.PCoA.Unifrac, type="samples", color = "Burning_intensity")+ 
  theme_bw()+
  ggtitle(label = "PCoA - Unweigthed Unifrac distance - ITS ")+
  geom_point(colour = "black", size = 4.5) +
  geom_point(size = 4)+
  geom_text_repel(aes(label=Sample_ID), color="black")+
  theme(axis.text=element_text(size=12))+ # size of lables, axis titles
  theme(axis.title.x = element_text(size=15), axis.title.y = element_text(size=15))+
  theme(legend.text=element_text(size=15),
        legend.title=element_text(size=15),
        legend.background = element_rect(size=0.5, linetype="solid", 
                                         colour ="black"))+
  scale_color_manual(values = c("#9cccd8","#ffb53f","#f26419","#931125"))+
  geom_line(aes(group = Plot))
PCoA_Unifrac_B



ggarrange(PCA_Ait, PCoA_Unifrac,
          ncol=2, nrow=1, 
          common.legend = TRUE, legend="bottom")


#### I WUnifrac distance ####

ord.PCoA.wUnifrac <- ordinate(ps.unmerged.tree, method="PCoA", distance = "wunifrac")

PCoA_wUnifrac <- plot_ordination(ps.unmerged.tree, ord.PCoA.wUnifrac, type="samples", color = "Plot",shape="Burning_intensity")+ 
  theme_bw()+
  ggtitle(label = "PCoA - Weigthed Unifrac distance - ITS ")+
  geom_point(colour = "black", size = 4.5) +
  geom_point(size = 4)+
  geom_text_repel(aes(label=Sample_ID), color="black")+
  theme(axis.text=element_text(size=12))+ # size of lables, axis titles
  theme(axis.title.x = element_text(size=15), axis.title.y = element_text(size=15))+
  theme(legend.text=element_text(size=15),
        legend.title=element_text(size=15),
        legend.background = element_rect(size=0.5, linetype="solid", 
                                         colour ="black"))+
  geom_line(aes(group = Plot))
PCoA_wUnifrac

PCoA_Unifrac_B <- plot_ordination(ps.unmerged.tree, ord.PCoA.Unifrac, type="samples", color = "Burning_intensity")+ 
  theme_bw()+
  ggtitle(label = "PCoA - Unweigthed Unifrac distance - ITS ")+
  geom_point(colour = "black", size = 4.5) +
  geom_point(size = 4)+
  geom_text_repel(aes(label=Sample_ID), color="black")+
  theme(axis.text=element_text(size=12))+ # size of lables, axis titles
  theme(axis.title.x = element_text(size=15), axis.title.y = element_text(size=15))+
  theme(legend.text=element_text(size=15),
        legend.title=element_text(size=15),
        legend.background = element_rect(size=0.5, linetype="solid", 
                                         colour ="black"))+
  scale_color_manual(values = c("#9cccd8","#ffb53f","#f26419","#931125"))+
  geom_line(aes(group = Plot))
PCoA_Unifrac_B





