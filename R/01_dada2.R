# This script follows the DADA2 pipeline tutorial to infer amplicon sequences from
# raw amplicon sequencing data.
# Original DADA2 tutorial available at https://benjjneb.github.io/dada2/tutorial.html

# Before running the script, download the raw amplicon sequence fastq files.
# Recommended: rename raw fastq files to remove sequencing prefixes.

# DADA2 (v. 1.12.1) pipeline

rm(list = ls())

library(dada2)
packageVersion("dada2")  # Should be v. 1.12.1

# Enter path to directory containing raw fastq files
path <- "PATH/DIRECTORY"
list.files(path)  # Lists raw fastq files

# Check number of samples
length(list.files(path))/2

fnFs <- sort(list.files(path, pattern = "_R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2.fastq", full.names = TRUE))

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Plot quality
plotQualityProfile(fnFs[1:10])  # Quality plots of forward reads in first 10 samples
plotQualityProfile(fnRs[1:10])  # Quality plots of reverse reads in first 10 samples

# Enter forward and reverse cut-off lengths based on above quality profile plots
forward_read_lengths <- INTEGER
reverse_read_lengths <- INTEGER

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Enter forward and reverse primer sequences in nucleotides
(forward_primer_length <- nchar("F_PRIMER_SEQUENCE"))  # Forward primer
(reverse_primer_length <- nchar("R_PRIMER_SEQUENCE"))  # Reverse primer

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft = c(forward_primer_length, reverse_primer_length), truncLen = c(forward_read_lengths, reverse_read_lengths), maxN = 0, maxEE = c(2,2), truncQ = 2, rm.phix = TRUE, compress = TRUE, multithread = TRUE)

#Learn error rates
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

# Sample inference
dadaFs <- dada(filtFs, err = errF, multithread = TRUE, pool = TRUE)
dadaRs <- dada(filtRs, err = errR, multithread = TRUE, pool = TRUE)
dadaFs[[1]]

# Merge paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)
head(mergers[[1]])

# Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)  # 118 x 21518
save(seqtab, file = "SEQTAB.rda")

# Save table summarizing fragment length frequencies
(nchar_table <- table(nchar(getSequences(seqtab))))
write.csv(nchar_table, file = "SEQ_LENGTH_FREQUENCIES_POOLED.csv", row.names = FALSE)

# Remove chimeras
seqtab_nochim_pooled <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE, verbose = TRUE)
dim(seqtab_nochim_pooled)
sum(seqtab_nochim_pooled)/sum(seqtab)

# Save seqtab.nochim as .rda file
save(seqtab_nochim_pooled, file = "SEQTAB_NOCHIM_POOLED.rda")

# Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab_nochim_pooled))
colnames(track) <- c("denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
write.csv(track, file = "DADA2_TRACK_READS_QAQC_POOLED.csv")

# Assign taxonomy
taxa_pooled <- assignTaxonomy(seqtab_nochim_pooled, "PATH/DB", taxLevels = c("kingdom","supergroup","division","class","order","family","genus","species"), multithread = TRUE)
save(taxa_pooled, file = "TAXA_POOLED.rda")
