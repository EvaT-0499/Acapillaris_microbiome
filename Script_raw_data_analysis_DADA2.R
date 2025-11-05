# Script for raw data treatment
# Eva Tanneau
# inspired from the DADA2 pipeline tutorial (https://benjjneb.github.io/dada2/tutorial.html)
# citation:
bibentry(bibtype="article", title="DADA2: High-resolution sample inference from Illumina amplicon data", author=c("Benjamin J Callahan", "Paul J McMurdie", "Michael J Rosen", "Andrew W Han", "Amy Jo A Johnson", "Susan P Holmes"), journal="Nature Methods", volume=13, pages="581-583", year=2016, doi="10.1038/nmeth.3869")

# Load packages
library(dada2); packageVersion("dada2")
library(ShortRead); packageVersion("ShortRead")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
library(phyloseq); packageVersion("phyloseq")
library(data.table); packageVersion("data.table")
library(stringr); packageVersion("stringr")
library(pander); packageVersion("pander")

# Example for ITS 2022 dataset
#Open path with raw data - Each library is first treatment separately
path.1 <- "/ITS_lib1_2022"
head(list.files(path.1))

fnFs.1 <- sort(list.files(path.1, pattern = "_R1_001.fastq.gz", full.names = TRUE))
fnRs.1 <- sort(list.files(path.1, pattern = "_R2_001.fastq.gz", full.names = TRUE))


path.2 <- "/ITS_lib2_2022"
head(list.files(path.2))

fnFs.2 <- sort(list.files(path.2, pattern = "_R1_001.fastq.gz", full.names = TRUE))
fnRs.2 <- sort(list.files(path.2, pattern = "_R2_001.fastq.gz", full.names = TRUE))


path.3 <- "/ITS_lib3_2022"
head(list.files(path.3))

fnFs.3 <- sort(list.files(path.3, pattern = "_R1_001.fastq.gz", full.names = TRUE))
fnRs.3 <- sort(list.files(path.3, pattern = "_R2_001.fastq.gz", full.names = TRUE))






# Create path for filtered data
fnFs.1.filtR <- file.path("/ITS_lib1_2022/filtR", 
                          basename(fnFs.1)) 
fnRs.1.filtR <- file.path("/ITS_lib1_2022/filtR", 
                          basename(fnRs.1))

fnFs.2.filtR <- file.path("/ITS_lib2_2022/filtR", 
                          basename(fnFs.2)) 
fnRs.2.filtR <- file.path("/ITS_lib2_2022/filtR", 
                          basename(fnRs.2))

fnFs.3.filtR <- file.path("/ITS_lib3_2022/filtR", 
                          basename(fnFs.3)) 
fnRs.3.filtR <- file.path("/ITS_lib3_2022/filtR", 
                          basename(fnRs.3))


# filtering and trimming step
filterAndTrim(fnFs.1, fnFs.1.filtR, fnRs.1, fnRs.1.filtR, matchIDs=TRUE, maxN=0, multithread = TRUE)

filterAndTrim(fnFs.2, fnFs.2.filtR, fnRs.2, fnRs.2.filtR, matchIDs=TRUE, maxN=0, multithread = TRUE)

filterAndTrim(fnFs.3, fnFs.3.filtR, fnRs.3, fnRs.3.filtR, matchIDs=TRUE, maxN=0, multithread = TRUE)



# Check quality profile of the raw data
#lib 1
plotQualityProfile(fnFs.1.filtR[1:10])
plotQualityProfile(fnRs.1.filtR[1:10])

plotQualityProfile(fnFs.1.filtR, aggregate = T)
plotQualityProfile(fnRs.1.filtR, aggregate = T)

#lib 2
plotQualityProfile(fnFs.2.filtR[1:10])
plotQualityProfile(fnRs.2.filtR[1:10])

plotQualityProfile(fnFs.2.filtR, aggregate = T)
plotQualityProfile(fnRs.2.filtR, aggregate = T)

#lib 3
plotQualityProfile(fnFs.3.filtR[1:10])
plotQualityProfile(fnRs.3.filtR[1:10])

plotQualityProfile(fnFs.3.filtR, aggregate = T)
plotQualityProfile(fnRs.3.filtR, aggregate = T)




######################### Remove primers
# Primer sequence
FWD <- "CTTGGTCATTTAGAGGAAGTAA"  ## change accordingly for ITS or 16S
REV <- "GCTGCGTTCTTCATCGATGC"  ## change accordingly for ITS or 16S

# To ensure we have the right primers, and the correct orientation of the primers on the reads, we will verify the presence and orientation of these primers in the data:
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

# Count how often primers appear in the forward and reverse reads - considering all possible primer orientations - for a subset of files only.
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}


#lib 1
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.1.filtR[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.1.filtR[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.1.filtR[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.1.filtR[[1]]))

#lib 2
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.2.filtR[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.2.filtR[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.2.filtR[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.2.filtR[[1]]))

#lib 3
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.3.filtR[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.3.filtR[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.3.filtR[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.3.filtR[[1]]))


# Remove primers using cutadapt
# call cutadapt from the emplacement on the computer
cutadapt <- "/cutadapt"
system2(cutadapt, args = "--version") # Run shell commands from R


# Create path for data after primer removal
#lib 1
path.1.cut <- file.path(path.1, "cutadapt")
if(!dir.exists(path.1.cut)) dir.create(path.1.cut)
fnFs.1.cut <- file.path(path.1.cut, basename(fnFs.1))
fnRs.1.cut <- file.path(path.1.cut, basename(fnRs.1))

#lib 2
path.2.cut <- file.path(path.2, "cutadapt")
if(!dir.exists(path.2.cut)) dir.create(path.2.cut)
fnFs.2.cut <- file.path(path.2.cut, basename(fnFs.2))
fnRs.2.cut <- file.path(path.2.cut, basename(fnRs.2))

#lib 3
path.3.cut <- file.path(path.3, "cutadapt")
if(!dir.exists(path.3.cut)) dir.create(path.3.cut)
fnFs.3.cut <- file.path(path.3.cut, basename(fnFs.3))
fnRs.3.cut <- file.path(path.3.cut, basename(fnRs.3))



FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)

R1.flags <- paste("-g", FWD, "-a", REV.RC) 
R2.flags <- paste("-G", REV, "-A", FWD.RC) 





# Run cutadapt:
#lib 1
for(i in seq_along(fnFs.1)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-j", 4, # multithreading
                             "-o", fnFs.1.cut[i], "-p", fnRs.1.cut[i], # output files
                             fnFs.1.filtR[i], fnRs.1.filtR[i])) # input files
}

#lib 2
for(i in seq_along(fnFs.2)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-j", 4, # multithreading
                             "-o", fnFs.2.cut[i], "-p", fnRs.2.cut[i], # output files
                             fnFs.2.filtR[i], fnRs.2.filtR[i])) # input files
}

#lib 3
for(i in seq_along(fnFs.3)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-j", 4, # multithreading
                             "-o", fnFs.3.cut[i], "-p", fnRs.3.cut[i], # output files
                             fnFs.3.filtR[i], fnRs.3.filtR[i])) # input files
}



# Sanity check before and after removing primers (counting)
# primer sequences lib1
# before
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.1.filtR[[11]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.1.filtR[[11]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.1.filtR[[11]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.1.filtR[[11]]))
# after
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.1.cut[[11]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.1.cut[[11]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.1.cut[[11]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.1.cut[[11]]))

# primer sequences lib2
# before
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.2.filtR[[11]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.2.filtR[[11]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.2.filtR[[11]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.2.filtR[[11]]))
# after
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.2.cut[[11]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.2.cut[[11]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.2.cut[[11]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.2.cut[[11]]))

# primer sequences lib3
# before
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.3.filtR[[11]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.3.filtR[[11]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.3.filtR[[11]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.3.filtR[[11]]))
# after
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.3.cut[[11]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.3.cut[[11]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.3.cut[[11]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.3.cut[[11]]))



###################### Separate forward and reverse sequences

# Separate reverse and forward sequences based on the pattern in the file name:
#lib 1
cutFs.1 <- sort(list.files(path.1.cut, pattern = "_R1_001.fastq.gz", full.names = TRUE))
cutRs.1 <- sort(list.files(path.1.cut, pattern = "_R2_001.fastq.gz", full.names = TRUE))

#lib 2
cutFs.2 <- sort(list.files(path.2.cut, pattern = "_R1_001.fastq.gz", full.names = TRUE))
cutRs.2 <- sort(list.files(path.2.cut, pattern = "_R2_001.fastq.gz", full.names = TRUE))

#lib 3
cutFs.3 <- sort(list.files(path.3.cut, pattern = "_R1_001.fastq.gz", full.names = TRUE))
cutRs.3 <- sort(list.files(path.3.cut, pattern = "_R2_001.fastq.gz", full.names = TRUE))



# Extract sample names, assuming file names have format:
# 20_61_CON_01_L_R2_L001_R1_001.fastq.gz
# We need them later to extract the files that "still exist" after quality filtering
#lib 1
get.sample.name <- function(fname) strsplit(basename(fname), "_L001")[[1]][1]
sample.names.1 <- unname(sapply(cutFs.1, get.sample.name))
head(sample.names.1)

#lib 2
get.sample.name <- function(fname) strsplit(basename(fname), "_L001")[[1]][1]
sample.names.2 <- unname(sapply(cutFs.2, get.sample.name))
head(sample.names.2)

#lib 3
get.sample.name <- function(fname) strsplit(basename(fname), "_L001")[[1]][1]
sample.names.3 <- unname(sapply(cutFs.3, get.sample.name))
head(sample.names.3)





########################### Filtering = Cut sequences to keep only the good quality part

# Create the path for separating forward and revers sequence files that pass the quality filtering:
#lib 1
filtFs.1 <- file.path(path.1.cut, "filtered", basename(cutFs.1))
filtRs.1 <- file.path(path.1.cut, "filtered", basename(cutRs.1))

#lib 2
filtFs.2 <- file.path(path.2.cut, "filtered", basename(cutFs.2))
filtRs.2 <- file.path(path.2.cut, "filtered", basename(cutRs.2))

#lib 3
filtFs.3 <- file.path(path.3.cut, "filtered", basename(cutFs.3))
filtRs.3 <- file.path(path.3.cut, "filtered", basename(cutRs.3))



# Cut the sequences removing the bad quality part (according to the quality plots earlier, cutting in maxLen and minLen)
#lib 1
out.1 <- filterAndTrim(cutFs.1, filtFs.1, cutRs.1, filtRs.1, maxN = 0, maxEE = c(2, 2), 
                       truncQ = 2, minLen = 100, maxLen = 400, rm.phix = TRUE, compress = TRUE, multithread = TRUE)

head(out.1)
#Save results in a csv file
write.csv(out.1, file="Reads_in-out_quality_filtering_lib1.csv", sep = ",")

#lib 2
out.2 <- filterAndTrim(cutFs.2, filtFs.2, cutRs.2, filtRs.2, maxN = 0, maxEE = c(2, 2), 
                       truncQ = 2, minLen = 100, maxLen = 400, rm.phix = TRUE, compress = TRUE, multithread = TRUE)

head(out.2)
#Save results in a csv file
write.csv(out.2, file="Reads_in-out_quality_filtering_lib2.csv", sep = ",")

#lib 3
out.3 <- filterAndTrim(cutFs.3, filtFs.3, cutRs.3, filtRs.3, maxN = 0, maxEE = c(2, 2), 
                       truncQ = 2, minLen = 100, maxLen = 400, rm.phix = TRUE, compress = TRUE, multithread = TRUE)

head(out.3)
#Save results in a csv file
write.csv(out.3, file="Reads_in-out_quality_filtering_lib3.csv", sep = ",")


# Select in "exists" only the filtered sequence files (forward):
exists.1 <- file.exists(filtFs.1)
exists.2 <- file.exists(filtFs.2)
exists.3 <- file.exists(filtFs.3)





########################################## Error treatment

# Error learning
errF.1 <- learnErrors(filtFs.1[exists.1], nbases= 1e10, randomize = TRUE, multithread = TRUE)
errR.1 <- learnErrors(filtRs.1[exists.1], nbases= 1e10, randomize = TRUE, multithread = TRUE)

errF.2 <- learnErrors(filtFs.2[exists.2], nbases= 1e10, randomize = TRUE, multithread = TRUE)
errR.2 <- learnErrors(filtRs.2[exists.2], nbases= 1e10, randomize = TRUE, multithread = TRUE)

errF.3 <- learnErrors(filtFs.3[exists.3], nbases= 1e10, randomize = TRUE, multithread = TRUE)
errR.3 <- learnErrors(filtRs.3[exists.3], nbases= 1e10, randomize = TRUE, multithread = TRUE)

# Save the error rate results:
saveRDS(errF.1, "/error_rate_F_lib1.rds")
saveRDS(errR.1, "/error_rate_R_lib1.rds")

saveRDS(errF.1, "/error_rate_F_lib2.rds")
saveRDS(errR.1, "/error_rate_R_lib2.rds")

saveRDS(errF.1, "/error_rate_F_lib3.rds")
saveRDS(errR.1, "/error_rate_R_lib3.rds")

# Plot error
plotErrors(errF.1, nominalQ=TRUE)
plotErrors(errR.1, nominalQ=TRUE)

plotErrors(errF.2, nominalQ=TRUE)
plotErrors(errR.2, nominalQ=TRUE)

plotErrors(errF.3, nominalQ=TRUE)
plotErrors(errR.3, nominalQ=TRUE)








####################################### Deprecation, merging and chimeras removal

# Dereplicate identical reads
derepFs.1 <- derepFastq(filtFs.1[exists.1], verbose=TRUE)
derepRs.1 <- derepFastq(filtRs.1[exists.1], verbose=TRUE)

derepFs.2 <- derepFastq(filtFs.2[exists.2], verbose=TRUE)
derepRs.2 <- derepFastq(filtRs.2[exists.2], verbose=TRUE)

derepFs.3 <- derepFastq(filtFs.3[exists.3], verbose=TRUE)
derepRs.3 <- derepFastq(filtRs.3[exists.3], verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs.1) <- sample.names.1[exists.1]
names(derepRs.1) <- sample.names.1[exists.1]

names(derepFs.2) <- sample.names.2[exists.2]
names(derepRs.2) <- sample.names.2[exists.2]

names(derepFs.3) <- sample.names.3[exists.3]
names(derepRs.3) <- sample.names.3[exists.3]

# ASV inference (denoising step):
dadaFs.1 <- dada(derepFs.1, err = errF.1, multithread = TRUE)
saveRDS(dadaFs.1, "/dada_F_lib1.rds")

dadaRs.1 <- dada(derepRs.1, err = errR.1, multithread = TRUE)
saveRDS(dadaRs.1, "/dada_R_lib1.rds")

# Sample inference lib 2
dadaFs.2 <- dada(derepFs.2, err = errF.2, multithread = TRUE)
saveRDS(dadaFs.2, "/dada_F_lib2.rds")

dadaRs.2 <- dada(derepRs.2, err = errR.2, multithread = TRUE)
saveRDS(dadaRs.2, "/dada_R_lib2.rds")

# Sample inference lib 3
dadaFs.3 <- dada(derepFs.3, err = errF.3, multithread = TRUE)
saveRDS(dadaFs.3, "/dada_F_lib3.rds")

dadaRs.3 <- dada(derepRs.3, err = errR.3, multithread = TRUE)
saveRDS(dadaRs.3, "/dada_R_lib3.rds")






# Merge forward and reverse reads
mergers.1 <- mergePairs(dadaFs.1, derepFs.1, dadaRs.1, derepRs.1, verbose=TRUE)
mergers.2 <- mergePairs(dadaFs.2, derepFs.2, dadaRs.2, derepRs.2, verbose=TRUE)
mergers.3 <- mergePairs(dadaFs.3, derepFs.3, dadaRs.3, derepRs.3, verbose=TRUE)

# Inspect the merger data.frame from the first sample of lib1
head(mergers.1[[1]])
head(mergers.2[[1]])
head(mergers.3[[1]])

# Save the merged sequences
saveRDS(mergers.1, "/merged-reads_lib1.rds")
saveRDS(mergers.2, "/merged-reads_lib2.rds")
saveRDS(mergers.3, "/merged-reads_lib3.rds")





# Generate sequence table from the merged reads
seqtab.1 <- makeSequenceTable(mergers.1)
dim(seqtab.1)

seqtab.2 <- makeSequenceTable(mergers.2)
dim(seqtab.2)

seqtab.3 <- makeSequenceTable(mergers.3)
dim(seqtab.3)

# Save sequence table
saveRDS(seqtab.1, "/seqtab.1.rds")
saveRDS(seqtab.2, "/seqtab.2.rds")
saveRDS(seqtab.3, "/seqtab.3.rds")





# Merge libraries
st.all <- mergeSequenceTables(seqtab.1, seqtab.2, seqtab.3)
saveRDS(st.all, "/seqtab.all.rds")

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(st.all, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(st.all)
(sum(st.all)-sum(seqtab.nochim))/(sum(st.all)/100)

#Save final sequence table as .rds file
saveRDS(seqtab.nochim, "/seqtab.nonchim.rds")




################################## Recap of the number of sequences at each step

# Get sequence length
length.dis <- data.frame(table(nchar(getSequences(seqtab.nochim))))

# Plot the sequence length
ggplot(length.dis, aes(Var1, Freq)) + geom_col() + ggtitle("Sequence length distribution") + xlab("length in bp") + ylab("frequency") + theme(axis.text.x = element_text(size=5, angle=90))


# Select unique sequences
getN <- function(x) sum(getUniques(x))

# Create table with summary of the remaining number of sequences after each step
# Save as csv file
#lib 1
out.1 <- read.table("Reads_in-out_quality_filtering_lib1.csv", header=T, sep = ",")
track.1 <- cbind(out.1, sapply(dadaFs.1, getN), sapply(dadaRs.1, getN), sapply(mergers.1, getN), rowSums(seqtab.1))
colnames(track.1) <- c("name", "input", "filtered", "denoisedF", "denoisedR", "merged", 
                       "nonchim")
head(track.1)
write.csv(track.1, file="Reads_quality_filtering_all-steps_lib1.csv", sep = ",", row.names = F)

#lib 2
out.2 <- read.table("Reads_in-out_quality_filtering_lib2.csv", header=T, sep = ",")
track.2 <- cbind(out.2, sapply(dadaFs.2, getN), sapply(dadaRs.2, getN), sapply(mergers.2, getN), rowSums(seqtab.2))
colnames(track.2) <- c("name", "input", "filtered", "denoisedF", "denoisedR", "merged", 
                       "nonchim")
head(track.2)
write.csv(track.2, file="Reads_quality_filtering_all-steps_lib2.csv", sep = ",", row.names = F)

#lib 3
out.3 <- read.table("Reads_in-out_quality_filtering_lib3.csv", header=T, sep = ",")
track.3 <- cbind(out.3, sapply(dadaFs.3, getN), sapply(dadaRs.3, getN), sapply(mergers.3, getN), rowSums(seqtab.3))
colnames(track.3) <- c("name", "input", "filtered", "denoisedF", "denoisedR", "merged", 
                       "nonchim")
head(track.3)
write.csv(track.3, file="Reads_quality_filtering_all-steps_lib03.csv", sep = ",", row.names = F)





############################### Creation of the physeq object (that will be used for the next analysis)

# load reference file for taxonomy (UNITE for ITS, Silva for 16S)
unite.ref <- "/sh_general_release_dynamic_10.05.2020.fasta"

# Assign taxonomy (according to the right database) to the sequence table
taxa <- assignTaxonomy(seqtab.nochim, unite.ref, multithread = TRUE, tryRC = TRUE)

taxa.print <- taxa  # Removing sequence row names for display only
rownames(taxa.print) <- NULL
head(taxa.print)

saveRDS(taxa, "/taxa_ITS_2022.rds")





# Get sample names of all samples that made it through the filters and turn into data.frame
samples.out <- rownames(seqtab.nochim)
samples.out.frame <- data.frame(samples.out)
colnames(samples.out.frame) <- c("Sample")

# Load "condition table" containing the sample details
samples.det <- read.csv("/condition_table.csv", sep = ";")

# merge samples.out and sample detail table - to only keep the details for the files that passed all filters
samdf <- merge(samples.out.frame,samples.det, by = "Sample")
rownames(samdf) <- samples.out
# turn 61 and 62 in Block from numeric into characters
samdf$Block <- as.character(samdf$Block)





# Construct phyloseq object (physeq or ps) assembling the taxonomy, the sequences and the condition table
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))

# Save phyloseq object as .rds file
saveRDS(ps, "/ps_raw_ITS_2022.rds")






# Add the sequencing depth in the physeq information
ps@sam_data$depth <- sample_sums(ps)

# plot sequencing depth before pruning
sdt = data.table(as(sample_data(ps), "data.frame"),
                 TotalReads = sample_sums(ps), keep.rownames = TRUE)
pSeqDepth = ggplot(sdt, aes(depth)) + geom_histogram(bins=40) + ggtitle("Sequencing Depth")
pSeqDepth




nsamples(ps) # total number of sample remaining
table(sample_data(ps)$Block) # number of sample in each block
table(sample_data(ps)$Treatment) #  number of sample per vegetation practices
table(sample_data(ps)$Material) # number of sample per plant tissue

# total number of "sequences"
sum(as.data.frame(sample_sums(ps)))
sum(as.data.frame(taxa_sums(ps)))




# Select only the samples with more than 100 reads remaining
ps <- prune_samples(sample_sums(ps) >= 100, ps) # 23 data sets removed

nsamples(ps) # total number of sample remaining
table(sample_data(ps)$Block) # number of sample in each block
table(sample_data(ps)$Treatment) #  number of sample per vegetation practices
table(sample_data(ps)$Material) # number of sample per plant tissue

# total number of "sequences"
sum(as.data.frame(sample_sums(ps)))
sum(as.data.frame(taxa_sums(ps)))

# plot sequencing depth
sdt = data.table(as(sample_data(ps), "data.frame"),
                 TotalReads = sample_sums(ps), keep.rownames = TRUE)
pSeqDepth = ggplot(sdt, aes(depth)) + geom_histogram(bins=40) + ggtitle("Sequencing Depth")
pSeqDepth




# Write taxa names = ASV sequences into dna
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)

# Write DNA sequences into refseq slot
ps <- merge_phyloseq(ps, dna)
# change taxa names to ASV plus running number
names.ps <- as.data.frame(cbind(taxa_names(ps),
                                paste0("ASV", seq(ntaxa(ps)))),
                          stringsAsFactors=F)
colnames(names.ps) <- c('seq', 'asv')
taxa_names(ps) <- names.ps$asv

# Save sequences as csv file
tab <- data.frame(refseq(ps))
write.csv(tab, file="/ITS_sequences_2022.csv", sep = ",", quote=F)

# Inspect refseq slot
ps
head(refseq(ps))




# Rarefaction of the samples according to the smallest size of samples
# This step transform your relative data in absolute data (comparable between samples)
set.seed(1); .Random.seed
ps_rarefied <- rarefy_even_depth(ps, sample.size=1*min(sample_sums(ps)))





