#### GETTING READY ####
R.version #4.1.1 
rm(list = ls())
gc()

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#                            DADA 2 script                                   #                                   
#                                                                            #
#    This script runs the full DADA2 bioinformatic pipeline using ITS data   #                  
#                                                                            #
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#

##### Packages #####
library(dada2); packageVersion("dada2") #1.20.0
library(ShortRead); packageVersion("ShortRead")#1.50.0
library(Biostrings); packageVersion("Biostrings")#2.60.2

setwd("/storage/Simon_morvan/Reads/ch2/ITS/02_DADA2_obj_mock")
path <- ("/storage/Simon_morvan/Reads/ch2/ITS/02_DADA2_obj_mock")  # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)


#### Identify primers  ####

# Forward and reverse fastq filenames have format: SAMPLENAME_R1.fastq and SAMPLENAME_R2.fastq
fnFs <- sort(list.files(path, pattern="MockFun_ITS_R1.fastq", full.names = TRUE)) # Fforwards reads
fnRs <- sort(list.files(path, pattern="MockFun_ITS_R2.fastq", full.names = TRUE)) # Reverse reads
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq


# Primers
FWD <- "GATGAAGAACGYAGYRA" # ITS3KYO2
REV <- "TCCTCCGCTTATTGATATGC" #ITS4

# Checking we have the correct primers 

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients


# Removing reads with N 
fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)


primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

cutadapt <- "/usr/bin/cutadapt" # CHANGE ME to the cutadapt path on your machine
system2(cutadapt, args = "--version") # Run shell commands from R
# 1.15


# Output filenames for the cutadapt-ed files
path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

# Running cutadapt
FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 
# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}

# Sanity check # 
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))
# Primers and there complements/reverse complements shouldn't be detected in your reads


# Forward and reverse cutadapt-ed fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "_R1.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_R2.fastq.gz", full.names = TRUE))

sample.names <- sapply(strsplit(basename(fnFs), "_R"), `[`, 1) # les noms d'échantillons sont les mêmes pour reverse et forwards seul le R1 / R2 chnge
save(sample.names,file="sample.names.RData")

## Filter and Trim 
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))

out <- filterAndTrim(cutFs, filtFs, 
                     cutRs, filtRs, 
                     maxN = 0, 
                     maxEE = c(2, 2), 
                     truncQ = 2, 
                     minLen = 50, # discard reads shorter than 50bp
                     rm.phix = TRUE, 
                     compress = TRUE, 
                     multithread = TRUE)  
save(out,file="out.RData")
head(out)






#### Learn the Error Rates ####
errF <- learnErrors(filtFs, 
                    multithread=TRUE,
                    randomize = TRUE,
                    verbose=TRUE) # For forward reads

save(errF,file="errF.RData")


errR <- learnErrors(filtRs, 
                    multithread=TRUE,
                    randomize = TRUE,
                    verbose=TRUE) # For reverse reads
save(errR,file="errR.RData")



#### Sample Inference ####

dadaFs <- dada(filtFs, 
               err=errF, 
               pool = "FALSE", 
               multithread=TRUE,
               verbose=2)
save(dadaFs,file="dadaFs.RData")
rm(list = ls(pattern="dadaFs")) 


dadaRs <- dada(filtRs, 
               err=errR, 
               pool = "FALSE", 
               multithread=TRUE,
               verbose=2)
save(dadaRs,file="dadaRs.RData")
rm(list = ls(pattern="dadaRs")) 

# Inspecting the dada-class object returned by dada:
# dadaFs[[1]]
# dadaRs[[x]]

# The DADA2 algorithm inferred x real sequence variants from the x unique sequences in the first sample. 
# There is much more to the dada-class return object than this (see help("dada-class") for some info), 
# including multiple diagnostics about the quality of each inferred sequence variant, 
# but that is beyond the scope of an introductory tutorial.

load("dadaFs.RData")
load("dadaRs.RData")

#### Merge paired reads ####
mergers <- mergePairs(dadaFs, filtFs, 
                      dadaRs, filtRs, 
                      minOverlap = 12,
                      verbose=TRUE)
save(mergers,file="mergers.RData")


mergers_with_reject<- mergePairs(dadaFs, filtFs, 
                                 dadaRs, filtRs, 
                                 minOverlap = 12,
                                 returnRejects = TRUE,
                                 verbose=TRUE)
save(mergers_with_reject,file="mergers_with_reject.RData")
# Contains the rejected sequences


# Inspect the merger data.frame from the first sample
head(mergers[[1]])


#### Construct sequence table ####

# We can now construct a sequence table of our mouse samples, 
# a higher-resolution version of the OTU table produced by traditional methods.
seqtab <- makeSequenceTable(mergers)
save(seqtab,file="seqtab.RData")


#### Remove chimeras ####
# The core dada method removes substitution and indel errors, but chimeras remain. 
# Fortunately, the accuracy of the sequences after denoising makes identifying chimeras simpler 
# than it is when dealing with fuzzy OTUs: all sequences which can be exactly reconstructed as a bimera
# (two-parent chimera) from more abundant sequences.

# Remove chimeric sequences:
?removeBimeraDenovo
seqtab.nochim <- removeBimeraDenovo(seqtab, 
                                    method="consensus",
                                    multithread=TRUE,
                                    verbose=TRUE)
save(seqtab.nochim,file="seqtab.nochim.RData")

seqtab.nochim.bin <- ifelse(seqtab.nochim>0,1,0) # transforms the number of occurrence of ASV in sample into presence/absence of ASV in sample
save(seqtab.nochim.bin,file="seqtab.nochim.bin.RData")


# Track
getN <- function(x) sum(getUniques(x))
track <- cbind(out, getN(dadaFs), getN(dadaRs), getN(mergers), rowSums(seqtab.nochim),rowSums(seqtab.nochim.bin))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim","ASV#")
rownames(track) <- sample.names
track <- as.data.frame(track)
track$ID <- row.names(track)
track$ID <- substr(track$ID,28,nchar(track$ID))
row.names(track) <- track$ID
head(track)
save(track,file="track.RData")



## Assign taxonomy 

unite_fun_ref <- "/storage/Simon_morvan/Reads/ch2/ITS/sh_general_release_dynamic_10.05.2021.fasta"  # CHANGE ME to location on your machine
taxa_unite_fun <- assignTaxonomy(seqtab.nochim, unite_fun_ref, multithread = TRUE, tryRC = TRUE)
setwd("/storage/Simon_morvan/Reads/ch2/ITS/02_DADA2_obj_mock")
save(taxa_unite_fun,file="Mock_taxa_unite_fun.RData")

unite_euk_ref <- "/storage/Simon_morvan/Reads/ch2/ITS/sh_general_release_dynamic_all_10.05.2021.fasta"
taxa_unite_euk<- assignTaxonomy(seqtab.nochim, unite_euk_ref, multithread = TRUE, tryRC = TRUE)
setwd("/storage/Simon_morvan/Reads/ch2/ITS/02_DADA2_obj_mock")
save(taxa_unite_euk,file="Mock_taxa_unite_euk.RData")



#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#                            THE END OF DADA2                                #                                   
#                                                                            #                                                                                 #
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#

library(dada2)
library(phyloseq)
setwd("/Users/simonmorvan/Desktop/UdeM/Projet/Ch2/ITS_Data_Analysis/02_DADA2_obj_mock/")

load("seqtab.nochim.RData")
ASV <- otu_table(t(seqtab.nochim), taxa_are_rows = TRUE)



load("Mock_taxa_unite_fun.RData")
for (i in 1:7){
  taxa_unite_fun[,i]<-substring(taxa_unite_fun[,i],4)
}
TAX <- tax_table(as.matrix(taxa_unite_fun))

ps.mock.ITS <- phyloseq(ASV,TAX)

Mock <- psmelt(ps.mock.ITS)

write.table(Mock,
            dec=".",
            sep=";",
            row.names= F,
            "/Users/simonmorvan/Desktop/UdeM/Projet/Ch2/ITS_Data_Analysis/02_DADA2_obj_mock/Mock.csv")
# CSV for downstream analysis of the mock community

