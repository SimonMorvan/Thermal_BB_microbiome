#### GETTING READY ####
R.version #4.1.1 
rm(list = ls())
gc()

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#                            DADA 2 script                                   #                                   
#                                                                            #
#    This script runs the full DADA2 bioinformatic pipeline using 16S mock   #                  
#                                                                            #
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#

# This script was run a the remote server using only the mock community fastq files.
# The objective was to prevent spurious ASV discovery in the mock
# which could be due to their presence in other samples 

##### Packages #####
library(dada2); packageVersion("dada2") #1.20.0
library(ShortRead); packageVersion("ShortRead")#1.50.0
library(Biostrings); packageVersion("Biostrings")#2.60.2

setwd("/storage/Simon_morvan/Reads/ch2/16S/02_DADA2_obj_mock")
path <- ("/storage/Simon_morvan/Reads/ch2/16S/02_DADA2_obj_mock")  # Change me to the directory containing the fastq files after unzipping.
list.files(path)


#### Identify primers  ####

# Forward and reverse fastq filenames have format: SAMPLENAME_R1.fastq and SAMPLENAME_R2.fastq
fnFs <- sort(list.files(path, pattern="MockBac_16S_R1.fastq", full.names = TRUE)) # Forwards reads
fnRs <- sort(list.files(path, pattern="MockBac_16S_R2.fastq", full.names = TRUE)) # Reverse reads
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq

# Primers
# Stagered primer 341F
FP1 <- "CCTACGGGNGGCWGCAG" 
FP2 <- "TCCTACGGGNGGCWGCAG"
FP3 <- "ACCCTACGGGNGGCWGCAG"
FP4 <- "CTACCTACGGGNGGCWGCAG"

# Stagered primer 805R
RP1 <- "GACTACHVGGGTATCTAATCC"
RP2 <- "TGACTACHVGGGTATCTAATCC"
RP3 <- "ACGACTACHVGGGTATCTAATCC"
RP4 <- "CTAGACTACHVGGGTATCTAATCC"


allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}

FP1.orients <- allOrients(FP1)
FP2.orients <- allOrients(FP2)
FP3.orients <- allOrients(FP3)
FP4.orients <- allOrients(FP4)

RP1.orients <- allOrients(RP1)
RP2.orients <- allOrients(RP2)
RP3.orients <- allOrients(RP3)
RP4.orients <- allOrients(RP4)

# Removing reads with N 
fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

# Check presence of every versions of primers 
rbind(FWD.ForwardReads = sapply(FP1.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FP1.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(RP1.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(RP1.orients, primerHits, fn = fnRs.filtN[[1]]))

rbind(FWD.ForwardReads = sapply(FP2.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FP2.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(RP2.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(RP2.orients, primerHits, fn = fnRs.filtN[[1]]))

rbind(FWD.ForwardReads = sapply(FP3.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FP3.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(RP3.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(RP3.orients, primerHits, fn = fnRs.filtN[[1]]))

rbind(FWD.ForwardReads = sapply(FP4.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FP4.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(RP4.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(RP4.orients, primerHits, fn = fnRs.filtN[[1]]))

#CUTADAPT
cutadapt <- "/usr/bin/cutadapt" # CHANGE ME to the cutadapt path on your machine
system2(cutadapt, args = "--version") # Run shell commands from R
# 1.15

# Running cutadapt
# Reverse complements 
FP1.RC <- dada2:::rc(FP1)
FP2.RC <- dada2:::rc(FP2)
FP3.RC <- dada2:::rc(FP3)
FP4.RC <- dada2:::rc(FP4)

RP1.RC <- dada2:::rc(RP1)
RP2.RC <- dada2:::rc(RP2)
RP3.RC <- dada2:::rc(RP3)
RP4.RC <- dada2:::rc(RP4)

# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FP1, "-g", FP2,"-g", FP3, "-g", FP4, "-a", RP1.RC,"-a", RP2.RC,"-a", RP3.RC,"-a", RP4.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", RP1, "-G", RP2,"-G", RP3,"-G", RP4,"-A", FP1.RC,"-A", FP2.RC,"-A", FP3.RC,"-A", FP4.RC) 


# Output filenames for the cutadapt-ed files
path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))


# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}


# Sanity check # 
# Primers and there complements/reverse complements shouldn't be detected in your reads
rbind(FWD.ForwardReads = sapply(FP1.orients, primerHits, fn = fnFs.cut[[1]]), 
      FWD.ReverseReads = sapply(FP1.orients, primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(RP1.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(RP1.orients, primerHits, fn = fnRs.cut[[1]]))

rbind(FWD.ForwardReads = sapply(FP2.orients, primerHits, fn = fnFs.cut[[1]]), 
      FWD.ReverseReads = sapply(FP2.orients, primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(RP2.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(RP2.orients, primerHits, fn = fnRs.cut[[1]]))

rbind(FWD.ForwardReads = sapply(FP3.orients, primerHits, fn = fnFs.cut[[1]]), 
      FWD.ReverseReads = sapply(FP3.orients, primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(RP3.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(RP3.orients, primerHits, fn = fnRs.cut[[1]]))

rbind(FWD.ForwardReads = sapply(FP4.orients, primerHits, fn = fnFs.cut[[1]]), 
      FWD.ReverseReads = sapply(FP4.orients, primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(RP4.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(RP4.orients, primerHits, fn = fnRs.cut[[1]]))

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
                     truncLen = c(270,240),
                     maxEE = c(2,2), 
                     truncQ = 2, 
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

seqtab.nochim <- removeBimeraDenovo(seqtab, 
                                    method="consensus",
                                    multithread=TRUE,
                                    verbose=TRUE)
save(seqtab.nochim,file="seqtab.nochim.RData")

seqtab.nochim.bin <- ifelse(seqtab.nochim>0,1,0) # transforms the number of occurence of ASV in sample into presence/absnece of ASV in sample
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

taxa <- assignTaxonomy(seqtab.nochim, "/storage/Simon_morvan/Reads/ch2/16S/silva_nr99_v138.1_wSpecies_train_set.fa.gz", multithread = TRUE, tryRC = TRUE)
save(taxa,file="16SMocktaxa.RData")



#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#                            THE END OF DADA2                                #                                   
#                                                                            #                                                                                 #
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#

# This part of the script was run on my personal computer 
# after transfering the files from the remote desktop

library(dada2)
library(phyloseq)
setwd("/Users/simonmorvan/Desktop/UdeM/Projet/Ch2/16S_Data_Analysis/02_DADA2_obj_mock/")

load("seqtab.nochim.RData")
ASV <- otu_table(t(seqtab.nochim), taxa_are_rows = TRUE)


load("16SMocktaxa.RData")
# for (i in 1:7){
#   taxa_unite_fun[,i]<-substring(taxa_unite_fun[,i],4)
# }
TAX <- tax_table(as.matrix(taxa))

ps.mock.16S <- phyloseq(ASV,TAX)

Mock <- psmelt(ps.mock.16S)

write.table(Mock,
            dec=".",
            sep=";",
            row.names= F,
            "/Users/simonmorvan/Desktop/UdeM/Projet/Ch2/16S_Data_Analysis/02_DADA2_obj_mock/Mock.csv")
# CSV file to inspect mock community