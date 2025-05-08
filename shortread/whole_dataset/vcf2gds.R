

# Libraries
library(data.table)
library(tidyverse)
library(SeqArray)
library(SNPRelate)
library(doParallel)

# Working directory
setwd("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/raw_vcf/")

# Executable in command line
out <-  c("combined.filtsnps10bpindels_snps_filthighmiss.qual.pass.ann.gds")
in.vcf <- c("combined.filtsnps10bpindels_snps_filthighmiss.qual.pass.ann.vcf.gz")

# Register cores
doParallel::registerDoParallel(cores = 15)

# Convert VCF to GDS
seqVCF2GDS(in.vcf, out, parallel = 15, verbose = TRUE)