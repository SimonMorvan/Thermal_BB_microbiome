#### GETTING READY ####
R.version #4.1.1 
rm(list = ls())
gc()

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#                            DADA 2 script                                   #                                   
#                                                                            #
#    This script runs the full DADA2 bioinformatic pipeline using 16S data   #                  
#                                                                            #
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#

# This script was run on a remote server. 

##### Packages #####
library(dada2); packageVersion("dada2") #1.20.0
library(ShortRead); packageVersion("ShortRead")#1.50.0
library(Biostrings); packageVersion("Biostrings")#2.60.2

setwd("/storage/Simon_morvan/Reads/ch2/16S/02_DADA2_obj")
path <- ("/storage/Simon_morvan/Reads/ch2/16S/")  # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)

#### Identify primers  ####

# Forward and reverse fastq filenames have format: SAMPLENAME_R1.fastq and SAMPLENAME_R2.fastq
fnFs <- sort(list.files(path, pattern="_R1.fastq", full.names = TRUE)) # Forward reads
fnRs <- sort(list.files(path, pattern="_R2.fastq", full.names = TRUE)) # Reverse reads
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

AGGATTAGATACCCBDGTAGTCTAG

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

# plotQualityProfile(filtFs[1:3])
# plotQualityProfile(filtRs[1:3])


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

# plotErrors(errF, nominalQ=TRUE))
# plotErrors(errR, nominalQ=TRUE))

#### Sample Inference ####

dadaFs <- dada(filtFs, 
               err=errF, 
               pool = "pseudo", 
               multithread=TRUE,
               verbose=2)
save(dadaFs,file="dadaFs.RData")
rm(list = ls(pattern="dadaFs")) 


dadaRs <- dada(filtRs, 
               err=errR, 
               pool = "pseudo", 
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
save(mergers_with_reject,file="mergers_with_reject.RData") # Shows which sequences have been rejected

#### Construct sequence table ####

# We can now construct a sequence table of our mouse samples, 
# a higher-resolution version of the OTU table produced by traditional methods.
seqtab <- makeSequenceTable(mergers)
save(seqtab,file="seqtab.RData")


# Remove chimeric sequences:

seqtab.nochim <- removeBimeraDenovo(seqtab, 
                                    method="pooled",
                                    multithread=TRUE,
                                    verbose=TRUE)
save(seqtab.nochim,file="seqtab.nochim.RData")


# Track
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
track <- as.data.frame(track)
track$ID <- row.names(track)
track$ID <- substr(track$ID,28,nchar(track$ID))
row.names(track) <- track$ID
head(track)
save(track,file="track_270_240_cutadapt.RData")

#### Taxonomy ####


taxa <- assignTaxonomy(seqtab.nochim, "/storage/Simon_morvan/Reads/ch2/16S/silva_nr99_v138.1_wSpecies_train_set.fa.gz", multithread = TRUE, tryRC = TRUE)
save(taxa,file="taxa.RData")

# END on remote server


