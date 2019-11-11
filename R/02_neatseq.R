# neatseq: Combine, analyze, and visualize amplicon sequence count and taxonomy data from DADA2
# Program created by Rebecca Garner

# Load libraries
library(tidyverse)
library(seqinr)

# Load seqtab.nochim (ASV table) output by DADA2
# **Enter file path and replace SEQTAB_NOCHIM.rda with ASV table R object name**
load("SEQTAB_NOCHIM.rda")

# Load taxa (taxonomy table) output by DADA2
# **Enter file path and replace TAXA.rda with taxonomy table R object name**
load("TAXA.rda")

# Row names are amplicon DNA sequences
# **Replace TAXA with loaded taxonomy table name**
taxonomy <- TAXA %>%
  as_tibble(rownames = "sequence")

# Number of unique ASVs
(n_asv <- n_distinct(taxonomy$sequence))

# Arbitrarily assign ASV codes to unique sequences
# **Consider adding a prefix to the asv code relating it to the amplicon sequencing dataset**
taxonomy_asv <- taxonomy %>%
  mutate(asv_code = paste0("ASV_", str_pad(string = 1:n_asv, width = nchar(n_asv), side = "left", pad = "0")))

# Reformat ASV table into long-format
# **Replace SEQTAB_NOCHIM with loaded ASV table name**
seq_melt <- SEQTAB_NOCHIM %>%
  as_tibble(rownames = "sample_id") %>%
  gather("sequence", "n_reads", -sample_id) %>%
  filter(n_reads > 0)

# Combine ASV and taxonomy tables
asv_melt <- seq_melt %>%
  left_join(taxonomy_asv, by = "sequence")

# Inspect singleton ASV sequences
# **Edit taxonomic ranks according to taxonomy assigned by DADA2 classifier
asv_melt %>%
  group_by(asv_code, kingdom, supergroup, division, class, order, family, genus, species) %>%
  summarize(n_reads_dataset = sum(n_reads)) %>%
  filter(n_reads_dataset == 1)

# # Generate .fasta file for taxa unassigned at supergroup rank
# unassigned_supergroup_seqs <- asv_melt %>%
#   filter(is.na(supergroup)) %>%
#   distinct(asv_code, sequence)
# 
# write.fasta(sequences = as.list(unassigned_supergroup_seqs$sequence),
#             names = unassigned_supergroup_seqs$asv_code,
#             file.out = "unassigned_supergroups.fasta")
