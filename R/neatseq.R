# neatseq: Combine, analyze, and visualize amplicon sequence count and taxonomy data from DADA2
# Program created by Rebecca Garner

# Load libraries
library(tidyverse)

# Set working directory to amplicon sequencing project directory
setwd("<WORKING_DIRECTORY>")

# Load seqtab.nochim (ASV table) output by DADA2
load("<SEQTAB_NOCHIM.rda>")

# Load taxa (taxonomy table) output by DADA2
load("<TAXA.rda>")
