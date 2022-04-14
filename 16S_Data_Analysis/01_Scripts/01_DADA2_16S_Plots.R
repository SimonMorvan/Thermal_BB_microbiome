#### GETTING READY ####
R.version #4.1.1 
rm(list = ls())
gc()

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#                            DADA 2 plots                                    #                                   
#                                                                            #
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#

# This script was run on my personal computer after transfering the files from the remote server

library(dada2); packageVersion("dada2") #1.20.0

setwd("/Users/simonmorvan/Desktop/UdeM/Projet/Ch2/16S_Data_Analysis/02_DADA2_obj/")
load("sample.names.RData")
filt_path <- ("/Users/simonmorvan/Desktop/UdeM/Projet/Ch2/16S_Data_Analysis/02_DADA2_obj/Filtered/")# Place filtered files in filtered subdirectory
list.files(filt_path)
filtFs <- file.path(filt_path,paste0(sample.names,"_F_filt.fastq.gz")) # Forwards reads
filtRs <- file.path(filt_path,paste0(sample.names,"_R_filt.fastq.gz")) # Reverse reads

plotQualityProfile(filtFs[1:32],aggregate = TRUE)
plotQualityProfile(filtRs[1:32],aggregate = TRUE)


load("errF.RData")
plotErrors(errF, nominalQ=TRUE)

load("errR.RData")
plotErrors(errR, nominalQ=TRUE)
