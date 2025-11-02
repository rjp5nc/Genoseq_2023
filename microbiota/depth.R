


#module load gcc/11.4.0  openmpi/4.1.4 icu R/4.3.1
R


library(vcfR)
library(adegenet)
library(ape)
library(dplyr)
library(proxy)
library(stringr)
library(ggplot2)


vcf <- read.vcfR("/scratch/rjp5nc/microbiota/chlorella/gvcf/chlorella_genotype.vcf.gz")

# Extract DP (Depth) values from the VCF file
dp_matrix <- extract.gt(vcf, element = "DP", as.numeric = TRUE)

# Calculate average depth per variant (row-wise mean)
avg_depth_per_sample <- colMeans(dp_matrix, na.rm = TRUE)

avg_depth_per_sample

write.csv(avg_depth_per_sample, "/scratch/rjp5nc/microbiota/chlorella/chlorella_depth.csv",
          row.names = FALSE)