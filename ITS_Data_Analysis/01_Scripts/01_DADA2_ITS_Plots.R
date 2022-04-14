#### GETTING READY ####
R.version #4.1.1 
rm(list = ls())
gc()

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#                            DADA 2 plots                                    #                                   
#                                                                            #
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#

# This script was run on my personal computer after transferring the files from the remote server

library(dada2); packageVersion("dada2") #1.20.0

load("sample.names.RData")
filt_path <- ("/Users/simonmorvan/Desktop/UdeM/Projet/Ch2/ITS_Data_Analysis/02_DADA2_obj/filtered/")# Place filtered files in filtered subdirectory
list.files(filt_path)
filtFs <- file.path(filt_path,paste0(sample.names,"_R1.fastq.gz")) # Forwards reads
filtRs <- file.path(filt_path,paste0(sample.names,"_R2.fastq.gz")) # Reverse reads

# Samples filtered QP
plotQualityProfile(filtFs[c((1:26),(28:33))],aggregate=TRUE)
plotQualityProfile(filtRs[c((1:26),(28:33))],aggregate=TRUE)

# Controls filtered QP
plotQualityProfile(filtFs[c(27,34)],aggregate=FALSE)
plotQualityProfile(filtRs[c(27,34)],aggregate=FALSE)


load("errF.RData")
plotErrors(errF, nominalQ=TRUE)

load("errR.RData")
plotErrors(errR, nominalQ=TRUE)
