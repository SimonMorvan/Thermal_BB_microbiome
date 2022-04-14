# Simon Morvan
# Sept 2021
R.version # 4.1.1

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#                            Data and visualization                          #                                   
#                                                                            #
#   This script aims to :                                                    # 
#                         - visualize DADA2 output (track object)            #                
#                         - to filter unwanted reads                        #                             
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#

Dir <- "/Users/simonmorvan/Desktop/UdeM/Projet/Ch2/16S_Data_Analysis/02_DADA2_obj/"
setwd(Dir)

#### Packages ####

library(phyloseq)
library(ggplot2)
library(ape)
library(dplyr)
library(Hmisc)
library(scales)
library(reshape2)
library(ggplot2)


#### Track reads through the pipeline ####
# As a final check of our progress, we’ll look at the number of reads that made it through 
# each step in the pipeline:

load("out.Rdata")
load("mergers.Rdata")
load("dadaFs.Rdata")
load("seqtab.Rdata")
load("seqtab.nochim.Rdata")
seqtab.nochim.bin <- ifelse(seqtab.nochim>0,1,0)# number of ASVs per sample
save(seqtab.nochim.bin,file="seqtab.nochim.bin.RData")
load("seqtab.nochim.bin.Rdata")
load("sample.names.Rdata")

getN <- function(x) sum(getUniques(x))
track <- cbind(out,  
               sapply(mergers, getN), 
               rowSums(seqtab.nochim),
               rowSums(seqtab.nochim.bin)) # number of ASVs per sample
colnames(track) <- c("Input", "Filtered","Merged","Non-chim","#ASVpersample")
track<-as.data.frame(track)


# Simplifying and ordering the sample names
track$ID <- rownames(track)
track$ID <- substr(track$ID,28,nchar(track$ID))
track$ID <- gsub("_16S", "", track$ID)
track$ID <- gsub("_R1.fastq.gz","",track$ID)

ordered_names <- c("1A","1B","2A","2B","3A","3B","4A","4B","5A","5B","6A","6B","7A","7B","8A","8B","9A","9B","10A","10B","11A","11B","12A","12B","13A","13B","14A","14B","15A","15B","16A","16B","BlankpcrCES","BLK1","MockBac")
track<- track[match(ordered_names, track$ID), ]  
row.names(track)<- track$ID 
save(track,file="track.RData")

# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)


###### Track visualization #######

load("track.RData") 

gtrack<- track[,c(1:4,6)]
colnames(gtrack) <- c("Filtering loss","Merging loss", "Chimera loss", "Retained reads","ID")
lgtrack <- melt(gtrack, id.vars="ID")
# long format

lgtrack$ID <- factor(lgtrack$ID, levels =c("1A","1B","2A","2B","3A","3B","4A","4B","5A","5B","6A","6B","7A","7B","8A","8B","9A","9B","10A","10B","11A","11B","12A","12B","13A","13B","14A","14B","15A","15B","16A","16B","BlankpcrCES","BLK1","MockBac"))

hist_track <- ggplot(lgtrack ,aes(x=ID, y=as.numeric(value), fill=variable)) +
  ggtitle("Track visualisation - 16S ") +
  geom_bar(stat="identity", position ="identity") + 
  theme_classic() + # Theme
  theme(plot.title = element_text(color="black", size=14, face="bold", hjust =0.5))+ # Title
  theme(panel.grid = element_blank(), panel.border = element_blank()) + # Removes the border and grid
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  # theme(legend.title = element_blank())+ #  delete legend title
  labs(fill = "Bio-informatic step")+
  scale_y_continuous(name="Abondance (reads)", expand = c(0,0), breaks=seq(from = 0, to = 200000, by = 10000)) 
# Changes axis title, deletes the space between the graph and the axis and sets the breaks. 
hist_track 



##### Track stats #####

load("track.RData") # Loads the track

track_stat <- track

lgtrack_samples <-melt(track_stat[1:32,], id.vars="ID") # Removes control (Blanks and Mock)
lgtrack_samples_stat <- lgtrack_samples %>%
                        group_by(variable)  %>%
                        summarise(min = min(value),
                                  max = max(value),
                                  mean = mean(value),
                                  sd = sd(value))

write.table(lgtrack_samples_stat,
            dec=".",
            sep=";",
            row.names= T,
            "/Users/simonmorvan/Desktop/UdeM/Projet/Ch2/16S_Data_Analysis/04_Figs&Analysis/Track_stats.csv")
# CSV table for data inspection

###______#####
#### Data filtering####

# We don't want to keep all of the ASVs identified by DADA2

# Taxonomy #
load(file="taxa.RData")
taxa[is.na(taxa)] <- "NA"
taxa <- as.data.frame(taxa)
table(taxa$Kingdom)
taxa.bac <- subset(taxa,taxa$Kingdom=="Bacteria") #Takes only the bac
table(taxa.bac$Kingdom)
# 12 ASV were assigned as Archeas 

table(taxa.bac$Order)
taxa.bac <- subset(taxa.bac,taxa.bac$Order!="Chloroplast")
# 22 Chloroplast

table(taxa.bac$Family)
taxa.bac <- subset(taxa.bac,taxa.bac$Family!="Mitochondria")
# 11 Mitochondria ASVs also removed 


# ASV.table # 
load("seqtab.nochim.RData")
ASV.table <- as.data.frame(t(seqtab.nochim))

# We'll also remove the singletons and doubletons 
ASV.nosingdoub <- subset(ASV.table,rowSums(ASV.table)>2)

ASV.sing <- subset(ASV.table,rowSums(ASV.table)==1) # 261 singletons
ASV.doub<- subset(ASV.table,rowSums(ASV.table)==2) # 615 doubletons

dim(ASV.table)[1]-dim(ASV.nosingdoub)[1]# 876 ASV were present only once or twice. 

# For plotting purposes we can remove the part of the sample names so we just have a number identifying the sample
colnames(ASV.nosingdoub)<- substr(colnames(ASV.nosingdoub),28,nchar(colnames(ASV.nosingdoub)))
colnames(ASV.nosingdoub) = gsub(pattern = "_R1.fastq.gz", replacement = "",colnames(ASV.nosingdoub))
colnames(ASV.nosingdoub) = gsub(pattern = "_16S", replacement = "",colnames(ASV.nosingdoub))


# Now we want to have matching sequence in both the Tax.table and ASV.table
# In order to do so we subset both tables so they have matching names. 

ASV.noSD <- subset(ASV.nosingdoub, rownames(ASV.nosingdoub)%in%rownames(taxa.bac))
# Deletes the non-bacterial ASVs from the sequence table (ASV.nosingdoub)

Tax.noSD <- subset(taxa.bac, rownames(taxa.bac)%in%rownames(ASV.noSD))
# Deletes singletons and doubletons ASVs from the taxonomy table (Tax.table)

dim(ASV.noSD)
dim(Tax.noSD)
# 4064 bacterial ASV with abundance > 2

identical(rownames(ASV.noSD),rownames(Tax.noSD))
# We have exactly the same ASVs in both the tax.table and the ASV.table

save(ASV.noSD,file="/Users/simonmorvan/Desktop/UdeM/Projet/Ch2/16S_Data_Analysis/03_RData/ASV.bac.noSD.RData")  
save(Tax.noSD,file="/Users/simonmorvan/Desktop/UdeM/Projet/Ch2/16S_Data_Analysis/03_RData/Tax.bac.noSD.RData")

