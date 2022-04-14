#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#                           DESeq 2 FUNGI                                 #                                   
#                                                                            #             
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
library(phyloseq)
library(DESeq2)

setwd("/Users/simonmorvan/Desktop/UdeM/Projet/Ch2/ITS_Data_Analysis/03_RData/")
list.files()

## USING ONLY NO BURN / HIGH BURN #### 
load("ps.merg.tree.RData")

ps_burn <- prune_samples(ps.merg.tree@sam_data$Burning_intensity %in% c("NoBurn","HighBurn"),ps.merg.tree)
ps_burn <- prune_taxa(taxa_sums(ps_burn)>0, ps_burn)
# 588 ASVs left 

### Convert Phyloseq to Deseq2 ### 
diagdds = phyloseq_to_deseq2(ps_burn, ~ Burning_intensity)
diagdds = estimateSizeFactors(diagdds, type = "poscounts")

diagdds = DESeq(diagdds,
                test = "Wald",
                fitType="local")
plotDispEsts(diagdds)

dispFit <-  diagdds@rowRanges@elementMetadata@listData$dispFit
dispGeneEst <- diagdds@rowRanges@elementMetadata@listData$dispGeneEst
test <- log(dispGeneEst)-log(dispFit)
med_abs_res <- median(abs(test),na.rm=TRUE)

diagdds2 = DESeq(diagdds,
                 test = "Wald",
                 fitType="parametric")
plotDispEsts(diagdds2)



dispFit2 <-  diagdds2@rowRanges@elementMetadata@listData$dispFit
dispGeneEst2 <- diagdds2@rowRanges@elementMetadata@listData$dispGeneEst
test2 <- log(dispGeneEst2)-log(dispFit2)
med_abs_res2 <- median(abs(test2),na.rm=TRUE)
# Take the fitType minimizing the median absolute residuals
# Here the parametric ! 

res = results(diagdds2, alpha=0.05) # results of the test FDR accepted = 5% 
plotMA(res)

res2 = res[order(res$padj, na.last=NA), ] # ordered and removes the ASVs which have padj = NAs 
plotMA(res2)

res <- as.data.frame(res)
summary(res)
summary(res2)

alpha = 0.05
sigtab = res2[(res2$padj < alpha), ] # Select significant padj
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps_burn)[rownames(sigtab), ], "matrix"))

otu <- as.data.frame(t(otu_table(ps_burn)))
# Taxa as rows, samples as columns

sigtab = cbind(sigtab,rowSums(otu[rownames(sigtab),]))
colnames(sigtab)[14] <- "Total_Abundance"

samp <- sample_data(ps_burn)
noBurn <- rownames(samp[samp$Burning_intensity == "NoBurn",])
highBurn <- rownames(samp[samp$Burning_intensity == "HighBurn",])

Ab_noBurn<- as.data.frame(rowSums(otu[,noBurn]))
Ab_highBurn <- as.data.frame(rowSums(otu[,highBurn]))

sigtab = cbind(sigtab,Ab_noBurn[rownames(sigtab),])
colnames(sigtab)[15] <- "Abundance_in_noBurn_samples"

sigtab = cbind(sigtab,Ab_highBurn[rownames(sigtab),])
colnames(sigtab)[16] <- "Abundance_in_highBurn_samples"






##### USING 4 GROUPS DIRECTLY ####

ASV_ab <- as.data.frame(otu_table(ps.merg.tree))
Metadata <- as.data.frame(sample_data(ps.merg.tree))


# set `1` as reference level
Metadata$Burning_intensity <- relevel(Metadata$Burning_intensity, ref = "NoBurn")

# run DESeq2

dds <- DESeqDataSetFromMatrix(
  countData = t(ASV_ab),
  colData = Metadata,
  design= ~ Burning_intensity)

dds = estimateSizeFactors(dds, type = "poscounts")

dds <- DESeq(dds, betaPrior = FALSE)

dds = DESeq(dds,
                test = "Wald",
                fitType="parametric")

# https://www.researchgate.net/post/How-to-interpret-results-of-DESeq2-with-more-than-two-experimental-groups
resultsNames(dds)
res <- results(dds, name="Burning_intensity_HighBurn_vs_NoBurn",alpha=0.05)
res2 = res[order(res$padj, na.last=NA), ] # ordered and removes the ASVs which have padj = NAs 

alpha = 0.05
sigtab2 = as.data.frame(res2[(res2$padj < alpha), ])# Select significant padj

sigtab2 = cbind(sigtab2, as(tax_table(ps.merg.tree)[rownames(sigtab2), ], "matrix"))

otu <- as.data.frame(t(otu_table(ps.merg.tree)))
# Taxa as rows, samples as columns

sigtab2 = cbind(sigtab2,rowSums(otu[rownames(sigtab2),]))
colnames(sigtab2)[14] <- "Total_Abundance"

samp <- sample_data(ps.merg.tree)
noBurn <- rownames(samp[samp$Burning_intensity == "NoBurn",])
highBurn <- rownames(samp[samp$Burning_intensity == "HighBurn",])

Ab_noBurn<- as.data.frame(rowSums(otu[,noBurn]))
Ab_highBurn <- as.data.frame(rowSums(otu[,highBurn]))

sigtab2 = cbind(sigtab2,Ab_noBurn[rownames(sigtab2),])
colnames(sigtab2)[15] <- "Abundance_in_noBurn_samples"

sigtab2 = cbind(sigtab2,Ab_highBurn[rownames(sigtab2),])
colnames(sigtab2)[16] <- "Abundance_in_highBurn_samples"

for (i in 1:dim(sigtab2)[1]){
  if (sigtab2$log2FoldChange[i]>0){
    sigtab2$regulation[i] = "More_abundant_in_Burn"}
  else {
    sigtab2$regulation[i] ="More_abundant_in_Ctrl"
  } 
}

summary(as.factor(sigtab2$regulation))



sigtabinf2 <- sigtab2
for (i in 1:dim(sigtabinf2)[1]){
  if (sigtabinf2$Abundance_in_noBurn_samples[i]==0){
    sigtabinf2$log2FoldChange[i]=25.5
  }
}

for (i in 1:dim(sigtabinf2)[1]){
  if (sigtabinf2$Abundance_in_highBurn_samples[i]==0){
    sigtabinf2$log2FoldChange[i]=-25.5
  }
}
library(ggplot2)
bp <- ggplot(sigtabinf2, aes(y=log2FoldChange, x=Genus, colour=regulation)) + 
  geom_jitter(height=0,width=0.3,shape=21,size=2) +
  scale_shape(solid=FALSE)+
  #geom_text_repel(data= sigtab, label=paste0(sigtab$Genus,"_",sigtab$Species),size=3)+
  geom_hline(yintercept = 0, colour="black",linetype=2)+
  geom_hline(yintercept = 24.7, colour= "red", linetype=2)+
  geom_hline(yintercept = -24.7, colour= "red", linetype=2)+
  theme(axis.title.x = element_blank(),
        axis.text.x= element_text(size=15, angle=90,hjust=1,vjust=0.5),
        axis.title.y = element_text(size=15),
        axis.text.y= element_text(size=15))+
  #theme(legend.position= "none")+
  scale_color_manual(values=c("steelblue1","firebrick"))

bp


######____#####

## TAX GLOM SP  ####
load("ps_sp_glom.RData")

## USING ONLY NO BURN / HIGH BURN #### 

ps_burn <- prune_samples(ps_sp_glom@sam_data$Burning_intensity %in% c("NoBurn","HighBurn"),ps_sp_glom)
ps_burn <- prune_taxa(taxa_sums(ps_burn)>0, ps_burn)
# 515 ASVs left 

t <- as.data.frame(tax_table(ps_sp_glom))

### Convert Phyloseq to Deseq2 ### 
diagdds = phyloseq_to_deseq2(ps_burn, ~ Burning_intensity)
diagdds = estimateSizeFactors(diagdds, type = "poscounts")
diagdds1 = DESeq(diagdds,
                test = "Wald",
                fitType="local",
                sfType="poscounts")

# plotDispEsts(diagdds)


dispFit1 <-  diagdds1@rowRanges@elementMetadata@listData$dispFit
dispGeneEst1 <- diagdds1@rowRanges@elementMetadata@listData$dispGeneEst
test1 <- log(dispGeneEst1)-log(dispFit1)
med_abs_res1 <- median(abs(test1),na.rm=TRUE)

diagdds2 = DESeq(diagdds,
                 test = "Wald",
                 fitType="parametric")
#plotDispEsts(diagdds2)

dispFit2 <-  diagdds2@rowRanges@elementMetadata@listData$dispFit
dispGeneEst2 <- diagdds2@rowRanges@elementMetadata@listData$dispGeneEst
test2 <- log(dispGeneEst2)-log(dispFit2)
med_abs_res2 <- median(abs(test2),na.rm=TRUE)

# Take the fitType minimizing the median absolute residuals
# Here the Parametric ! 

res = results(diagdds2, alpha=0.05) # results of the test FDR accepted = 5% 
res2 = res[order(res$padj, na.last=NA), ] # ordered and removes the ASVs which have padj = NAs 
plotMA(res2)
summary(res2)

alpha = 0.05
sigtab = res2[(res2$padj < alpha), ] # Select significant padj
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps_burn)[rownames(sigtab), ], "matrix"))


otu <- as.data.frame(t(otu_table(ps_burn)))
# Taxa as rows, samples as columns

sigtab = cbind(sigtab,rowSums(otu[rownames(sigtab),]))
colnames(sigtab)[14] <- "Total_Abundance"

samp <- sample_data(ps_burn)
noBurn <- rownames(samp[samp$Burning_intensity == "NoBurn",])
highBurn <- rownames(samp[samp$Burning_intensity == "HighBurn",])

Ab_noBurn<- as.data.frame(rowSums(otu[,noBurn]))
Ab_highBurn <- as.data.frame(rowSums(otu[,highBurn]))

sigtab = cbind(sigtab,Ab_noBurn[rownames(sigtab),])
colnames(sigtab)[15] <- "Abundance_in_noBurn_samples"

sigtab = cbind(sigtab,Ab_highBurn[rownames(sigtab),])
colnames(sigtab)[16] <- "Abundance_in_highBurn_samples"

# 3 ASVs identified 



##### USING 4 GROUPS DIRECTLY ####

ASV_ab <- as.data.frame(otu_table(ps_sp_glom))
Metadata <- as.data.frame(sample_data(ps_sp_glom))


# set "NoBurn" as reference level
Metadata$Burning_intensity <- relevel(Metadata$Burning_intensity, ref = "NoBurn")

# run DESeq2

dds <- DESeqDataSetFromMatrix(
  countData = t(ASV_ab),
  colData = Metadata,
  design= ~ Burning_intensity)

dds = estimateSizeFactors(dds, type = "poscounts")

dds1 = DESeq(dds,
            test = "Wald",
            fitType="parametric",
            sfType = "poscounts")


dispFit3 <-  dds1@rowRanges@elementMetadata@listData$dispFit
dispGeneEst3 <- dds1@rowRanges@elementMetadata@listData$dispGeneEst
test3 <- log(dispGeneEst3)-log(dispFit3)
med_abs_res3 <- median(abs(test3),na.rm=TRUE)


dds2 = DESeq(dds,
            test = "Wald",
            fitType="local",
            sfType = "poscounts")


dispFit4 <-  dds2@rowRanges@elementMetadata@listData$dispFit
dispGeneEst4 <- dds2@rowRanges@elementMetadata@listData$dispGeneEst
test4 <- log(dispGeneEst4)-log(dispFit4)
med_abs_res4 <- median(abs(test4),na.rm=TRUE)

# Here parametric lower as well

# https://www.researchgate.net/post/How-to-interpret-results-of-DESeq2-with-more-than-two-experimental-groups
resultsNames(dds1)
res <- results(dds1, name="Burning_intensity_HighBurn_vs_NoBurn",alpha=0.05)
res2 = res[order(res$padj, na.last=NA), ] # ordered and removes the ASVs which have padj = NAs 

alpha = 0.05
sigtab2 = as.data.frame(res2[(res2$padj < alpha), ])# Select significant padj

sigtab2 = cbind(sigtab2, as(tax_table(ps_sp_glom)[rownames(sigtab2), ], "matrix"))

otu <- as.data.frame(t(otu_table(ps_sp_glom)))
# Taxa as rows, samples as columns

sigtab2 = cbind(sigtab2,rowSums(otu[rownames(sigtab2),]))
colnames(sigtab2)[14] <- "Total_Abundance"

samp <- sample_data(ps_sp_glom)
noBurn <- rownames(samp[samp$Burning_intensity == "NoBurn",])
highBurn <- rownames(samp[samp$Burning_intensity == "HighBurn",])

Ab_noBurn<- as.data.frame(rowSums(otu[,noBurn]))
Ab_noBurn_samples <- as.data.frame(otu[,noBurn])

Ab_highBurn <- as.data.frame(rowSums(otu[,highBurn]))
Ab_highBurn_samples <- as.data.frame(otu[,highBurn])

sigtab2 <-  cbind(sigtab2,Ab_noBurn[rownames(sigtab2),])
colnames(sigtab2)[15] <- "Abundance_in_noBurn_samples"
sigtab2 <-  cbind(sigtab2,Ab_noBurn_samples[rownames(sigtab2),])

sigtab2 = cbind(sigtab2,Ab_highBurn[rownames(sigtab2),])
colnames(sigtab2)[20] <- "Abundance_in_highBurn_samples"
sigtab2 <-  cbind(sigtab2,Ab_highBurn_samples[rownames(sigtab2),])

for (i in 1:dim(sigtab2)[1]){
  if (sigtab2$log2FoldChange[i]>0){
    sigtab2$regulation[i] = "More_abundant_in_Burn"}
  else {
    sigtab2$regulation[i] ="More_abundant_in_Ctrl"
  } 
}

summary(as.factor(sigtab2$regulation))

# 12 ASVs found
# More_abundant_in_Burn More_abundant_in_Ctrl 
#       5                     7 

# Filtering 
filtering <- sigtab2[,c("1","14","6","9","2","8","12","15")]
filtering_prev <- as.data.frame(ifelse(filtering>0,1,0))
filtered_prev <- subset(filtering_prev,rowSums(filtering_prev[,1:4])>=3 | rowSums(filtering_prev[,5:8])>=3 )


sigtab3 <- subset(sigtab2,row.names(sigtab2)%in%row.names(filtered_prev))



sigtabinf2 <- sigtab3
for (i in 1:dim(sigtabinf2)[1]){
  if (sigtabinf2$Abundance_in_noBurn_samples[i]==0){
    sigtabinf2$log2FoldChange[i]=25.5
  }
}

for (i in 1:dim(sigtabinf2)[1]){
  if (sigtabinf2$Abundance_in_highBurn_samples[i]==0){
    sigtabinf2$log2FoldChange[i]=-25.5
  }
}


library(ggplot2)
bp <- ggplot(sigtabinf2, aes(x=log2FoldChange, y=paste(Genus,"sp."), fill=regulation)) + 
  geom_jitter(height=0.2,width=0,shape=21,size=3) +
  scale_shape(solid=FALSE)+
  geom_vline(xintercept = 0, colour="black",linetype=2)+
  geom_vline(xintercept = 24, colour= "red", linetype=2)+
  geom_vline(xintercept = -24, colour= "red", linetype=2)+
  labs(y="Genus", fill="Differentially abundant ASVs" )+
  theme(axis.title.x = element_text(size=15),
        axis.text.x= element_text(size=15, angle=90,hjust=1,vjust=0.5),
        axis.title.y = element_text(size=15),
        axis.text.y= element_text(size=15,face = "italic"),
        legend.text=element_text(size=15),
        legend.title=element_text(size=15),
        legend.position = "bottom",
        legend.direction = "vertical",
        legend.background = element_rect(size=0.5, linetype="solid", 
                                         colour ="black"))+
  #theme(legend.position= "none")+
  scale_fill_manual(labels = c("High intensity burn", "No burn"), values=c("#931125","#9cccd8"))

bp

