R.version # 4.1.1.

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#                               Phyloseq                                                                       
#                                                                            
#   This script aims to : - create phyloseq object
#                         - inspect controls
#                         - plot histograms based on taxonomy
#                         - plot rarefaction curves 
#                         - show NA per taxonomic rank                       
#                         - merge pseudoreplicates                                                   
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#


Dir <- "/Users/simonmorvan/Desktop/UdeM/Projet/Ch2/16S_Data_Analysis/"
setwd(Dir)

#### Packages ####
library(phyloseq)
library(ggplot2)
library(plyr)
library(dplyr)


#### Phyloseq object ######
list.files()
# Sequences
load(file = "03_RData/ASV.bac.noSD.RData")
ASV = otu_table(ASV.noSD, taxa_are_rows = TRUE)

# Taxonomy
load(file = "03_RData/Tax.bac.noSD.RData")
Taxo <- as.matrix(Tax.noSD[,1:7])
TAX = tax_table(Taxo)


# Phyloseq object
ps.bac = phyloseq(ASV, TAX)
ps.bac
sample_names(ps.bac)


# Splitting phyloseq object into two 

# Phyloseq object for controls
controls <- c("BlankpcrCES","BLK1","MockBac")
ps.bac.controls <- prune_samples((sample_names(ps.bac) %in% controls),ps.bac)
ps.bac.controls <- prune_taxa(taxa_sums(ps.bac.controls)>0, ps.bac.controls)
ASV_controls <- as.data.frame(otu_table(ps.bac.controls))
ASV_controls_bin <- ifelse(ASV_controls>0,1,0) 
colSums(ASV_controls_bin)
#    MockBac        BLK1 BlankpcrCES 
#       43          18           1 


# Common ASVs between samples & extraction blank 
ASV_exblank <- ASV_controls[2] # Extraction blank data
colnames(ASV_exblank)
ASV_exblank <- subset(ASV_exblank,rowSums(ASV_exblank)>0)
colSums(ASV_exblank)

ASV_samples <- ASV.noSD[ , -which(names(ASV.noSD) %in% controls)]
ASV_samples <- subset(ASV_samples,rowSums(ASV_samples)>0) # Delete ASvs not present in samples

common_ex_ASV <- intersect(rownames(ASV_samples), rownames(ASV_exblank))# Common ASVs between sample data set and exblank
# 7 ASV found in common in the samples dataset and the extraction blank

ASV_exblank_cont <- subset(ASV_exblank,rownames(ASV_exblank)%in%common_ex_ASV)
ASV_exblank_cont$BLK1
# 4    3    4 2053    7    1    5


ASV_samples_excont <- subset(ASV_samples,rownames(ASV_samples)%in%common_ex_ASV)
# ASV abundance = x in sample dataset present in x sample :

  Tax_exblank_cont<- subset(Tax.noSD,rownames(Tax.noSD)%in%common_ex_ASV) # Delete ASvs not present in samples
# ASV belongs to x taxa. 


# Common ASVs between samples & PCR blank 
ASV_pcrblank <- ASV_controls[3] # PCR blank data
colnames(ASV_pcrblank)
ASV_pcrblank <- subset(ASV_pcrblank,rowSums(ASV_pcrblank)>0)

ASV_samples <- ASV.noSD[ , -which(names(ASV.noSD) %in% controls)]
ASV_samples <- subset(ASV_samples,rowSums(ASV_samples)>0) # Delete ASvs not present in samples

common_pcr_ASV <- intersect(rownames(ASV_samples), rownames(ASV_pcrblank))# Common ASVs between sample data set and pcrblank
# 0 ASV found in both the samples dataset and the extraction blank


####_____#####
##### Controls taxonomy #####
Tax_controls <- psmelt(ps.bac.controls)
HowManyOrders <- dim(table(Tax_controls$order))

# https://medialab.github.io/iwanthue/
Palette <- c("#90682c",
             "#6963cd",
             "#86b93b",
             "#af60d3",
             "#52c05a",
             "#c549a8",
             "#4c8b38",
             "#d84681",
             "#62c394",
             "#d93b42",
             "#48b6d2",
             "#cd5f25",
             "#5176ba",
             "#b9ae39",
             "#9e9ae2",
             "#dc9636",
             "#bd76b4",
             "#697329",
             "#9c4463",
             "#36845f",
             "#e7705e",
             "#acae67",
             "#a74435",
             "#da986a",
             "#df828a")

#libr",ary(randomcoloR)
#palette <- distinctColorPalette(HowManyOrders)


#library(RColorBrewer)
#getPalette = colorRampPalette(brewer.pal(8, "Accent")) #8 being the number of colours in Set2
#OrderPalette = getPalette(HowManyOrders)

ggplot(Tax_controls, aes(x = Sample, y = Abundance, fill = order)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values=Palette)+
  theme(legend.position="bottom")+
  theme(axis.title.x = element_blank())

##### Sample taxonomy#####
# Samples
controls <- c("BlankpcrCES","BLK1","MockBac")
ps.bac.samples<- prune_samples(!(sample_names(ps.bac) %in% controls),ps.bac)
ps.bac.samples<- prune_taxa(!(taxa_names(ps.bac.samples) %in% rownames(ASV_exblank_cont)),ps.bac.samples) # Remove possible contamination
ps.bac.samples <- prune_taxa(taxa_sums(ps.bac.samples)>0,ps.bac.samples)

sample_names(ps.bac.samples)
sample.data <-read.csv2("/Users/simonmorvan/Desktop/UdeM/Projet/Ch2/16S_Data_Analysis/03_Metadata/metadata.csv",row.names = 1)
sample.data$Sample_ID<- as.factor(sample.data$Sample_ID)
sample.data$Plot <- as.factor(sample.data$Plot)
sample.data$Block <- as.factor(sample.data$Block)
sample.data$Burning_intensity <- factor(sample.data$Burning_intensity, levels = c("NoBurn", "LowBurn", "MedBurn","HighBurn"))


SDAT <- sample_data(sample.data)

ps.bac.samples <- merge_phyloseq(ps.bac.samples,SDAT)
save(ps.bac.samples,file="/Users/simonmorvan/Desktop/UdeM/Projet/Ch2/16S_Data_Analysis/03_RData/ps_bac_unmerged.RData")
# Phyloseq object containing only the samples but where the pseudoreplicates haven't been merged yet


Tax_samples <- psmelt(ps.bac.samples)
HowManyphyla <- dim(table(Tax_samples$Phylum))
Palette <- c("#d05286","#5ec36d","#a76ed1","#98b342","#6473db","#cd9c2e","#543585","#a0bb6a","#c26fc3","#42772b","#d471b2","#47bb8a","#cc453d","#33d4d1","#d74f67","#5e87d3","#cc6f37","#872861","#d3a257","#a03a49","#856f24","#dd7a67","#ccef47")



ggplot(Tax_samples, aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  facet_grid(.~Plot,drop=TRUE,scales="free_x")+
  scale_fill_manual(values=Palette)+
  theme(legend.position="bottom")+
  theme(axis.title.x = element_blank())

# Relab
ps.bac.samples.relab <- transform_sample_counts(ps.bac.samples,function(x) x/sum(x))# Relative abondance

Tax_samples_relab <- psmelt(ps.bac.samples.relab)
Palette <- c("#d05286","#5ec36d","#a76ed1","#98b342","#6473db","#cd9c2e","#543585","#a0bb6a","#c26fc3","#42772b","#d471b2","#47bb8a","#cc453d","#33d4d1","#d74f67","#5e87d3","#cc6f37","#872861","#d3a257","#a03a49","#856f24","#dd7a67","#ccef47")

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

setwd("/Users/simonmorvan/Desktop/UdeM/Projet/Ch2/16S_Data_Analysis/03_RData/")
load(file = "ps_bac_unmerged.RData")
ASV.ab <- as.data.frame(t(otu_table(ps.bac.samples)))

S <- specnumber(ASV.ab) # Observed number of ASVs per sample
max(S)
mean(S)
median(S)
min(S)

(raremax <- min(rowSums(ASV.ab))) # Minimum number of sequences (retrieves the abundance in the site that has the least sequences ) 
Srare <- rarefy(ASV.ab, raremax) # Gives the expected species richness in random subsamples of a definit size (raremax) 


Rarecurve_16S <- rarecurve(ASV.ab, 
                           step = 50, 
                           label = TRUE, 
                           sample = raremax, col ="blue", cex =0.5, 
                           xlab="Number of sequences", ylab="Number of ASV", 
                           cex.axis=1.3, cex.lab=1.4, cex.main=1.5)




####_____#####

##### NA Per taxonomy rank #### 
# This piece of code is counting the number of NA per taxonomic rank

load("ps_bac_unmerged.RData")
Tax <- as.data.frame(tax_table(ps.bac.samples))

mat=matrix(length(Tax[,1]),4,7)
# Creates a matrix, with 3 rows and 7 columns (Taxranks)
# Each of the cell countains the total amount of bacterial ASV noSD
mat<-as.data.frame(mat)
# Transform it into dataframe
colnames(mat)[1:7]<- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
rownames(mat)[1:4]<- c("ASV taxonomicly assigned ","Percentage (%)"," NA per rank","NA (sum)")
# specify rownames 
for (c in 1:7){
  for (i in 1:length(Tax[,1])){
    if ((Tax[i,c])=="NA")
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
save(NA.taxrank,file="NA.per.taxrank.RData")

library(reshape2)
plot.NA.taxrank <- NA.taxrank[c(1,4),]
rownames(plot.NA.taxrank)[1] <- "Taxonomy assigned"
rownames(plot.NA.taxrank)[2] <- "NA"
plot.NA.taxrank$ID <- rownames(plot.NA.taxrank)
lg.NA.tax <- melt(plot.NA.taxrank, id.vars="ID")


NA_hist <- ggplot(lg.NA.tax ,aes(x=variable, y=as.numeric(value), fill=ID)) +
  ggtitle("Depth of taxonomic assignment") +
  geom_bar(stat="identity", position ="stack") + 
  theme_classic() + # Theme
  theme(plot.title = element_text(color="black", size=14, face="bold", hjust =0.5))+ # Title
  theme(panel.grid = element_blank(), panel.border = element_blank()) + # Removes the border and grid
  theme(axis.ticks.length=unit(0.1,"cm")) + # Ticks size
  labs(fill = "Taxonomic status")+
  scale_x_discrete(name ="Taxonomic Rank", expand = c(0,0)) + 
  scale_y_continuous(name="ASVs", expand = c(0,0), breaks=seq(from = 0, to = 11000, by = 1000)) 
NA_hist 




#### _______ #####
#### Merging Pseudoreplicates ####
load(file ="03_RData/ps_bac_unmerged.RData")
test <- psmelt(ps.bac.samples)
sample_data(ps.bac.samples)[["Plot"]] <- as.character(sample_data(ps.bac.samples)[["Plot"]])
sample_data(ps.bac.samples)[["Burning_intensity"]]

ps.bac.merg <- merge_samples(ps.bac.samples, "Plot")

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
sample_data(ps.bac.merg)$Burning_intensity <- as.factor(sample_data(ps.bac.merg)$Burning_intensity)
sample_data(ps.bac.merg)$Burning_intensity <- revalue(sample_data(ps.bac.merg)$Burning_intensity, c("1"="NoBurn", "2" = "LowBurn", "3" = "MedBurn", "4" = "HighBurn"))
sample_data(ps.bac.merg)$Burning_intensity <- factor(sample_data(ps.bac.merg)$Burning_intensity, levels=c("NoBurn","LowBurn","MedBurn","HighBurn"))

sample_data(ps.bac.merg)<- subset(sample_data(ps.bac.merg), select=-c(Sample_ID))
sample_data(ps.bac.merg)$Plot <- as.factor(sample_data(ps.bac.merg)$Plot)
sample_data(ps.bac.merg)$Block <- as.factor(sample_data(ps.bac.merg)$Block)


save(ps.bac.merg,file="/Users/simonmorvan/Desktop/UdeM/Projet/Ch2/16S_Data_Analysis/03_RData/ps_bac_merg.RData")


##### Merged pseudoreplicate taxonomy #####
load("03_RData/ps_bac_merg.RData")
Tax_samples <- psmelt(ps.bac.merg)
HowManyphyla <- dim(table(Tax_samples$Phylum))
Palette <- c("#d05286","#5ec36d","#a76ed1","#98b342","#6473db","#cd9c2e","#543585","#a0bb6a","#c26fc3","#42772b","#d471b2","#47bb8a","#cc453d","#33d4d1","#d74f67","#5e87d3","#cc6f37","#872861","#d3a257","#a03a49","#856f24","#dd7a67","#ccef47")

ggplot(Tax_samples, aes(x = reorder(Sample, (as.numeric(as.character(Sample)))), y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values=Palette)+
  theme(legend.position="bottom")+
  theme(axis.title.x = element_blank())


# RelAb
ps.bac.merg.relab <- transform_sample_counts(ps.bac.merg,function(x) x/sum(x))# Relative abondance
Tax_samples_relab <- psmelt(ps.bac.merg.relab)

ggplot(Tax_samples_relab, aes(x = reorder(Sample, (as.numeric(as.character(Sample)))), y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  facet_grid(.~Burning_intensity,drop=TRUE,scales="free_x")+
  scale_fill_manual(values=Palette)+
  theme(legend.position="bottom")+
  theme(axis.title.x = element_blank())



