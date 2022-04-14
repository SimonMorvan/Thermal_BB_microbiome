# Simon Morvan
# Sept 2021
R.version # 4.1.1

Dir <- "/Users/simonmorvan/Desktop/UdeM/Projet/Ch2/ITS_Data_Analysis/02_DADA2_obj/"
setwd(Dir)

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#                            Data and visualization                          #                                   
#                                                                            #
#   This script aims to :                                                    # 
#                         - visualize DADA2 output (track object)            #                
#                         - to filter unwanted reads                                                                           #
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#

#### Packages####

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

# Simplifying and ordering the sample names
load("track.RData")
track$ID <- gsub("_ITS", "", track$ID)
ordered_names <- c("1A","1B","2A","2B","3A","3B","4A","4B","5A","5B","6A","6B","7A","7B","8A","8B","9A","9B","10A","10B","11A","11B","12A","12B","13A","13B","14A","14B","15A","15B","16A","16B","BlankpcrCES","BLK1")
track<- track[match(ordered_names, track$ID), ]  
row.names(track)<- track$ID 

save(track,file="track.RData")


gtrack<- track[,c(1,2,5,6,8)]
colnames(gtrack)[1] <- "Filtering loss"
colnames(gtrack)[2] <- "Merging loss"
colnames(gtrack)[3] <- "Chimera loss"
colnames(gtrack)[4] <- "Reads kept"
gtrack$ID <- factor(gtrack$ID,levels = c("1A","1B","2A","2B","3A","3B","4A","4B","5A","5B","6A","6B","7A","7B","8A","8B","9A","9B","10A","10B","11A","11B","12A","12B","13A","13B","14A","14B","15A","15B","16A","16B","BlankpcrCES","BLK1"))

lgtrack <- melt(gtrack, id.vars="ID")
# long format

###### Track visualization #######

hist_track <- ggplot(lgtrack ,aes(x=ID, y=as.numeric(value), fill=variable)) +
  ggtitle("Track visualisation - ITS ") +
  geom_bar(stat="identity", position ="identity") + 
  theme_classic() + # Theme
  theme(plot.title = element_text(color="black", size=14, face="bold", hjust =0.5))+ # Title
  theme(panel.grid = element_blank(), panel.border = element_blank()) + # Removes the border and grid
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        axis.title.x = element_blank())+
  # theme(legend.title = element_blank())+ #  delete legend title
  labs(fill = "Bio-informatic step")+
  scale_y_continuous(name="Abundance (reads)", expand = c(0,0), breaks=seq(from = 0, to = 200000, by = 10000)) 
# Changes axis title, deletes the space between the graph and the axis and sets the breaks. 
hist_track 



##### Track stats #####


load("track.RData") # Loads the track

track_stat <- track[,c(1,2,5,6:8)]

lgtrack_samples <-melt(track_stat[!row.names(track_stat)%in% c("BlankpcrCES","BLK1"),], id.vars="ID")

test <- subset(lgtrack_samples,lgtrack_samples$variable=="ASV#")
mean(test$value)

lgtrack_samples_stat <- lgtrack_samples %>%
                        group_by(variable)  %>%
                        summarise(min = min(value),
                                  max = max(value),
                                  sum = sum(value),
                                  mean = mean(value),
                                  sd = sd(value))

write.table(lgtrack_samples_stat,
            dec=".",
            sep=";",
            row.names= F,
            "/Users/simonmorvan/Desktop/UdeM/Projet/Ch2/ITS_Data_Analysis/04_Figs&Analysis/Track_stats.csv")
# CSV table for data inspection


###______#####
#### Data filtering####
# We don't want to keep all of the ASVs identified by DADA2

# Taxonomy #
load(file="taxa_unite_fun.RData") # Taxa obtained with fungal ref
load(file="taxa_unite_euk.RData") # Taxa obtained with eukaryots ref

for (i in 1:7){
  taxa_unite_euk[,i]<-substring(taxa_unite_euk[,i],4) # Removes the first 3 characters (k--, p-- etc)
  taxa_unite_fun[,i]<-substring(taxa_unite_fun[,i],4)
  }

taxa_unite_euk <- as.data.frame(taxa_unite_euk)
taxa_unite_fun <- as.data.frame(taxa_unite_fun)

taxa_unite_euk[] <- lapply(taxa_unite_euk, factor)
summary(taxa_unite_euk) # 784 Fungi on 913 obs. 

taxa_unite_fun[] <- lapply(taxa_unite_fun, factor)
summary(taxa_unite_fun) # 106 NAs as Phylum on 930 obs. 

incertain_Fungi <- subset(taxa_unite_fun,is.na(taxa_unite_fun$Phylum))
# 106  obs. obtained with the fungal ref don't have a Phylum
# Since the Fungal ref doesn't seem to have an out group, we want to check if these ASVs are indeed fungi. 

# Therefore we look what is there tax assignment when using the eukaryotic ref. database
coresp_euk <- subset(taxa_unite_euk, rownames(taxa_unite_euk)%in%rownames(incertain_Fungi))
coresp_euk[] <- lapply( coresp_euk, factor) 
summary(coresp_euk) # Only 1 ASVs is confirmed to be a Fungi, 
# most of the 106 ASVs are labelled as plants with the euk. ref. db. 


# Doing the analysis the other way around -> 
# checking which ASVs don't match as Fungi in the tax obtained with the euk. ref
incertain_Fungi2 <- subset(taxa_unite_euk,(taxa_unite_euk$Kingdom!="Fungi" | is.na(taxa_unite_euk$Kingdom)))
# 129 obs from Euk ref are classified different from Fungi (other Kingdom or NAs)
coresp_fun <-  subset(taxa_unite_fun, rownames(taxa_unite_fun)%in%rownames(incertain_Fungi2))
coresp_fun[] <- lapply( coresp_fun, factor) 
summary(coresp_fun)
# In tax from fungal ref, 105 of these 129 have NA in Phylum.
# Other 24 have a Phylum in taxa obtained with fungal ref. 

# Let's check if the 105 NA's found in coresp_fun 
# correspond to the 106 ASVs that did not match to a Fungal Kingdom 
# in corresp_euk
coresp_funNA <- subset(coresp_fun,is.na(coresp_fun$Phylum))
common_nonfungal <- intersect(rownames(coresp_funNA),rownames(coresp_euk))
diff <- subset(coresp_euk,!rownames(coresp_euk)%in%common_nonfungal)

# It is the case, therefore we'll remove the 105 observations that did not have a fungal Phylum 
# in the tax assigned by fungal database and which matched to Kingdoms different from fungi in the euk ref db. 
# We'll keep the taxonomy obtained from the fungal database. 
taxa.fun <- subset(taxa_unite_fun,!rownames(taxa_unite_fun)%in%rownames(coresp_funNA))



# ASV.table # 
load("seqtab.nochim.RData")
ASV.table <- as.data.frame(t(seqtab.nochim))

# We'll also remove the singletons and doubletons 
ASV.nosingdoub <- subset(ASV.table,rowSums(ASV.table)>2)

ASV.sing <- subset(ASV.table,rowSums(ASV.table)==1) # 0 singletons
ASV.doub<- subset(ASV.table,rowSums(ASV.table)==2) # 45 doubletons

dim(ASV.table)[1]-dim(ASV.nosingdoub)[1]# 45 ASV were present only once or twice. 

# For plotting purposes we can remove the SB_ so we just have the number identifying the sample
colnames(ASV.nosingdoub)<- substr(colnames(ASV.nosingdoub),28,nchar(colnames(ASV.nosingdoub)))
colnames(ASV.nosingdoub) = gsub(pattern = "_R1.fastq.gz", replacement = "",colnames(ASV.nosingdoub))# Removes SB_ from names
colnames(ASV.nosingdoub) = gsub(pattern = "_ITS", replacement = "",colnames(ASV.nosingdoub))# Removes SB_ from names



# Now we want to have matching sequence in both the Tax.table and ASV.table
# In order to do so we subset both tables so they have matching names. 

ASV.noSD <- subset(ASV.nosingdoub, rownames(ASV.nosingdoub)%in%rownames(taxa.fun))
# Deletes the non-fungal ASV from the sequence table (ASV.nosingdoub)

Tax.noSD <- subset(taxa.fun, rownames(taxa.fun)%in%rownames(ASV.noSD))
# Deletes ASVs that singletons and doubletons from the taxonomy table (Tax.table)

dim(ASV.noSD)
dim(Tax.noSD)
# 772 fungal ASV with abundance > 2

identical(rownames(ASV.noSD),rownames(Tax.noSD))
# We have exactly the same ASVs in both the tax.table and the ASV.table

asv_seqs <- rownames(ASV.noSD)
asv_headers <- vector(dim(ASV.noSD)[1], mode="character")
 for (i in 1:dim(ASV.noSD)[1]) {
   asv_headers[i] <- paste("ASV", i, sep="_")
 }
 
# making and writing out a fasta of our final ASV seqs:
asv_fasta <- cbind(asv_headers, asv_seqs)
write.table(asv_fasta,
            dec=".",
            sep=";",
            row.names= F,
            "/Users/simonmorvan/Desktop/UdeM/Projet/Ch2/ITS_Data_Analysis/03_RData/ASV.fasta.csv")



ASV.noSD.ab <- ASV.noSD
ASV.noSD.tax <- Tax.noSD
# giving our seq headers more manageable names (ASV_1, ASV_2...)
# row.names(ASV.noSD.ab) <- sub(">", "", asv_headers)
# rownames(ASV.noSD.tax) <- gsub(pattern=">", replacement="", x=asv_headers)



save(ASV.noSD.ab,file="/Users/simonmorvan/Desktop/UdeM/Projet/Ch2/ITS_Data_Analysis/03_RData/ASV.noSD.ab.RData")  
save(ASV.noSD.tax,file="/Users/simonmorvan/Desktop/UdeM/Projet/Ch2/ITS_Data_Analysis/03_RData/ASV.noSD.tax.RData")
