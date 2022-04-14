# Simon Morvan
# Sept 2O21
R.version # 4.1.1.
Dir <- "/Users/simonmorvan/Desktop/UdeM/Projet/Ch2/ITS_Data_Analysis/03_RData/"
setwd(Dir)
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#                               Phyloseq                                     #                                   
#                                                                            #
#   This script aims a diversity and ordinations                             #
#                                                                            #
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#

#### Packages ####
library(phyloseq)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(plyr)
library(dplyr)


#### Phyloseq object ######
list.files()
# Sequences
load(file = "ASV.noSD.ab.RData")
ASV = otu_table(ASV.noSD.ab, taxa_are_rows = TRUE)

# Taxonomy
load(file = "ASV.noSD.tax.RData")
Taxo <- as.matrix(ASV.noSD.tax[,1:7])
TAX = tax_table(Taxo)


# Phyloseq object
ps.fun = phyloseq(ASV, TAX)
ps.fun
sample_names(ps.fun)


# Splitting phyloseq object into two 


# Controls
controls <- c("BlankpcrCES","BLK1")
ps.fun.controls <- prune_samples((sample_names(ps.fun) %in% controls),ps.fun)
ps.fun.controls <- prune_taxa(taxa_sums(ps.fun.controls)>0, ps.fun.controls)
ASV_controls <- as.data.frame(otu_table(ps.fun.controls))
colSums(ASV_controls)
ASV_controls_bin <- ifelse(ASV_controls>0,1,0) 
colSums(ASV_controls_bin)


# Common ASVs between samples & extraction blank 
ASV_exblank <- ASV_controls["BLK1"] # Extraction blank data
ASV_exblank <- subset(ASV_exblank,rowSums(ASV_exblank)>0)


ASV_samples <- ASV.noSD.ab[ , -which(names(ASV.noSD.ab) %in% controls)]
ASV_samples <- subset(ASV_samples,rowSums(ASV_samples)>0) # Delete ASvs not present in samples

common_ex_ASV <- intersect(rownames(ASV_samples), rownames(ASV_exblank))# Common ASVs between sample data set and exblank
# 1 ASV found in common in the samples dataset and the extraction blank

ASV_exblank_cont <- subset(ASV_exblank,rownames(ASV_exblank)%in%common_ex_ASV)
# ASV abundance = 1 read in extraction blank

ASV_samples_excont <- subset(ASV_samples,rownames(ASV_samples)%in%common_ex_ASV)
rowSums(ASV_samples_excont)
ASV_samples_excont_bin <- ifelse(ASV_samples_excont>0,1,0)
rowSums(ASV_samples_excont_bin)
# ASV abundance = 437 in sample dataset present in 21 samples :

Tax_exblank_cont<- subset(ASV.noSD.tax,rownames(ASV.noSD.tax)%in%common_ex_ASV) 
# Troposporella monospora

####_____#####
##### Controls taxonomy#####
Tax_controls <- psmelt(ps.fun.controls)
HowManyOrders <- dim(table(Tax_controls$Order))

# https://medialab.github.io/iwanthue/
Palette <- c("#be73ab",
             "#5ab84e",
             "#cd51b0",
             "#a4b249",
             "#8161ca",
             "#d89b48",
             "#688bcd",
             "#cc5537",
             "#56c5a9",
             "#c8566d",
             "#2f8a72",
             "#8c6e32",
             "#53833f")

#libr",ary(randomcoloR)
#palette <- distinctColorPalette(HowManyOrders)


#library(RColorBrewer)
#getPalette = colorRampPalette(brewer.pal(8, "Accent")) #8 being the number of colours in Set2
#OrderPalette = getPalette(HowManyOrders)

ggplot(Tax_controls, aes(x = Sample, y = Abundance, fill = Order)) +
  geom_bar(stat = "identity",color="black") +
  #geom_bar(stat = "identity", position = "stack", color = "black")+
  scale_fill_manual(values=Palette)+
  theme(legend.position="bottom")+
  # guides(fill=guide_legend(nrow=))
  theme(axis.title.x = element_blank())


#### Mock taxonomy ### 
sample_names(ps.fun)
ps.mock <- prune_samples((sample_names(ps.fun) == "MockFun") ,ps.fun)
ps.mock <- prune_taxa(taxa_sums(ps.mock)>0, ps.mock)

Mock_taxa <- psmelt(ps.mock)
write.table(Mock_taxa,
            dec=".",
            sep=";",
            row.names= T,
            "/Users/simonmorvan/Desktop/UdeM/Projet/Ch2/ITS_Data_Analysis/03_RData/MockBac.csv")

  
  
  
  
##### Sample taxonomy#####
# Samples
controls <- c("BlankpcrCES","BLK1")
ps.fun.samples<- prune_samples(!(sample_names(ps.fun) %in% controls),ps.fun)# Remove samples 
ps.fun.samples<- prune_taxa(!(taxa_names(ps.fun.samples) %in% row.names(ASV_controls)),ps.fun.samples)
ps.fun.samples <- prune_taxa(taxa_sums(ps.fun.samples)>0,ps.fun.samples)



sample.data <-read.csv2("/Users/simonmorvan/Desktop/UdeM/Projet/Ch2/ITS_Data_Analysis/03_Metadata/metadata.csv",row.names = 1)
sample.data$Sample_ID<- as.factor(sample.data$Sample_ID)
sample.data$Plot <- as.factor(sample.data$Plot)
sample.data$Block <- as.factor(sample.data$Block)
sample.data$Burning_intensity <- factor(sample.data$Burning_intensity, levels = c("NoBurn", "LowBurn", "MedBurn","HighBurn"))

SDAT <- sample_data(sample.data)

ps.fun.samples <- merge_phyloseq(ps.fun.samples,SDAT)
save(ps.fun.samples,file="ps_fun_unmerged.RData")


Tax_samples <- psmelt(ps.fun.samples)
HowManyphyla <- dim(table(Tax_samples$Phylum))
Palette <- c("#6588cd",
             "#7ea342",
             "#a361c7",
             "#4aab83",
             "#d14066",
             "#c18e40",
             "#c26594")


ggplot(Tax_samples, aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  facet_grid(.~Plot,drop=TRUE,scales="free_x")+
  scale_fill_manual(values=Palette)+
  theme(legend.position="bottom")+
  theme(axis.title.x = element_blank())

# Relab
ps.fun.samples.relab <- transform_sample_counts(ps.fun.samples,function(x) x/sum(x))# Relative abondance

Tax_samples_relab <- psmelt(ps.fun.samples.relab)

## Comparing pseudoreplicates 
ggplot(Tax_samples_relab, aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  facet_grid(.~Plot,drop=TRUE,scales="free_x")+
  scale_fill_manual(values=Palette)+
  scale_y_continuous(name ="Relative abundance", labels = scales::percent, expand = c(0,0))+
  theme(legend.position="bottom")+
  theme(axis.title.x = element_blank())


## Effect of burning
ggplot(Tax_samples_relab, aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  facet_grid(.~Burning_intensity,drop=TRUE,scales="free_x")+
  scale_fill_manual(values=Palette)+
  scale_y_continuous(name ="Relative abundance", labels = scales::percent, expand = c(0,0))+
  theme(legend.position="bottom")+
  theme(axis.title.x = element_blank())


####_____#####
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#####                            Rarefaction curves                       ####                                   
#                                                                            #
#  This script allows to verify that the sequencing allows to show the full  #                  
#  microbial richness in each samples. Did we manage to capture the entire   #
#                 microbial diversity for every sample?                      #
#                                                                            #
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#

library(vegan)
library(ggplot2)

setwd("/Users/simonmorvan/Desktop/UdeM/Projet/Ch2/ITS_Data_Analysis/03_RData/")
load(file = "ps_fun_unmerged.RData")
ASV.ab <- as.data.frame(t(otu_table(ps.fun.samples)))
# Abundance table for fungi excluding singletons and doubletons 
# Samples as rows, taxa as columns 

S <- specnumber(ASV.ab) # Observed number of ASVs per sample
max(S)
mean(S)
median(S)
min(S) #BLK PCR 

(raremax <- min(rowSums(ASV.ab))) # Minimum number of sequences (retrieves the abundance in the site that has the least sequences ) 
Srare <- rarefy(ASV.ab, raremax) # Gives the expected species richness in random subsamples of a definit size (raremax) 


# RARECURVE #
Rarecurve_ITS <- rarecurve(ASV.ab, 
                           step = 50, 
                           label = TRUE, 
                           col ="blue", cex =0.5, 
                           xlab="Number of sequences", ylab="Number of ASVs", 
                           cex.axis=1.3, cex.lab=1.4, cex.main=1.5)



####_____#####
##### NA Per taxonomy rank #### 
# This piece of code is counting the number of NA per taxonomic rank

load("ps_fun_unmerged.RData")
Tax <- as.data.frame(tax_table(ps.fun.samples))

mat=matrix(length(Tax[,1]),4,7)
# Creates a matrix, with 3 rows and 7 columns (Taxranks)
# Each of the cell countains the total amount of bacterial ASV noSD
mat<-as.data.frame(mat)
# Transform it into dataframe
colnames(mat)[1:7]<- c("Domain","Phylum","Class","Order","Family","Genus","Species")
rownames(mat)[1:4]<- c("ASV taxonomicly assigned ","Percentage (%)"," NA per rank","NA (sum)")
# specify rownames 
for (c in 1:7){
  for (i in 1:length(Tax[,1])){
    if (is.na(Tax[i,c]))
      mat[,c] <- mat[,c] - 1
    as.numeric(mat[,c])
  }
}
mat[2,] <- round(((mat[1,]/length(Tax[,1]))*100),1)
mat[3,1]=0
mat[4,1]=0
for (i in 2:7){
  mat[3,i] =   mat[1,(i-1)] - mat[1,i]
  mat[4,i] = sum(mat[3,i],mat[4,i-1])
} 
mat <- as.data.frame(mat)
NA.taxrank <- mat
save(NA.taxrank,file="/Users/simonmorvan/Desktop/UdeM/Projet/Ch2/ITS_Data_Analysis/03_RData/NA.per.taxrank.RData")

library(reshape2)
plot.NA.taxrank <- NA.taxrank[c(1,4),]
rownames(plot.NA.taxrank)[1] <- "Taxonomy assigned"
rownames(plot.NA.taxrank)[2] <- "NA"
plot.NA.taxrank$ID <- rownames(plot.NA.taxrank)
lg.NA.tax <- melt(plot.NA.taxrank, id.vars="ID")


NA_hist <- ggplot(lg.NA.tax ,aes(x=variable, y=as.numeric(value), fill=ID)) +
  ggtitle("Depth of taxonomic assignment - ITS Data") +
  geom_bar(stat="identity", position ="stack") + 
  theme_classic() + # Theme
  theme(plot.title = element_text(color="black", size=14, face="bold", hjust =0.5))+ # Title
  theme(panel.grid = element_blank(), panel.border = element_blank()) + # Removes the border and grid
  theme(axis.ticks.length=unit(0.1,"cm")) + # Ticks size
  labs(fill = "Taxonomic status")+
  scale_x_discrete(name ="Taxonomic Rank", expand = c(0,0)) + 
  scale_y_continuous(name="ASVs", expand = c(0,0), breaks=seq(from = 0, to = 1000, by = 100)) 
NA_hist 


#### _______ #####
#### Merging Pseudoreplicates ####
load(file ="ps_fun_unmerged.RData")

sample_data(ps.fun.samples)[["Plot"]] <- as.character(sample_data(ps.fun.samples)[["Plot"]])
sample_data(ps.fun.samples)[["Burning_intensity"]]

ps.fun.merg <- merge_samples(ps.fun.samples, "Plot")

#Check
# OTUnames10 = names(sort(taxa_sums(ps.bac), TRUE)[1:10])#10 most abundant asv
# ps.bac10  = prune_taxa(OTUnames10,  ps.bac)# 
# Merged10 = prune_taxa(OTUnames10, ps.bac.merg)
# check <- as.data.frame((otu_table(ps.bac10)))
# check2 <- as.data.frame(t(otu_table(Merged10)))
# A <- rowSums(check)                       
# B <- rowSums(check2)                       
# identical(A,B)
# 
# check3 <- as.data.frame(otu_table(ps.bac))
# check4 <- as.data.frame(t(otu_table(ps.bac.merg)))


# Resetting Burning_intensity as factors and ordering them
sample_data(ps.fun.merg)$Burning_intensity <- as.factor(sample_data(ps.fun.merg)$Burning_intensity)
sample_data(ps.fun.merg)$Burning_intensity <- revalue(sample_data(ps.fun.merg)$Burning_intensity, c("1"="HighBurn", "2" = "LowBurn", "3" = "MedBurn", "4" = "NoBurn"))
sample_data(ps.fun.merg)$Burning_intensity <- factor(sample_data(ps.fun.merg)$Burning_intensity, levels=c("NoBurn","LowBurn","MedBurn","HighBurn"))


sample_data(ps.fun.merg)<- subset(sample_data(ps.fun.merg), select=-c(Sample_ID))
sample_data(ps.fun.merg)$Plot <- as.factor(sample_data(ps.fun.merg)$Plot)
sample_data(ps.fun.merg)$Block <- as.factor(sample_data(ps.fun.merg)$Block)


save(ps.fun.merg,file="ps.merg.RData")



##### Merged pseudoreplicate taxonomy #####
load("03_RData/ps.merg.RData")
Tax_samples <- psmelt(ps.fun.merg)
HowManyphyla <- dim(table(Tax_samples$Phylum))

Palette <- c("#6588cd",
             "#7ea342",
             "#a361c7",
             "#4aab83",
             "#d14066",
             "#c18e40",
             "#c26594")

Raw_ab<- ggplot(Tax_samples, aes(x = reorder(Sample, (as.numeric(as.character(Sample)))), y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values=Palette)+
  theme(legend.position="bottom")+
  theme(axis.title.x = element_blank())
Raw_ab


# RelAb
ps.fun.merg.relab <- transform_sample_counts(ps.fun.merg,function(x) x/sum(x))# Relative abondance
Tax_samples_relab <- psmelt(ps.fun.merg.relab)

ggplot(Tax_samples_relab, aes(x = reorder(Sample, (as.numeric(as.character(Sample)))), y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  facet_grid(.~Burning_intensity,drop=TRUE,scales="free_x")+
  scale_fill_manual(values=Palette)+
  theme(legend.position="bottom")+
  theme(axis.title.x = element_blank())

