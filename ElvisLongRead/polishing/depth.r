#module load gcc/11.4.0  openmpi/4.1.4 icu R/4.3.1

library(vcfR)
library(adegenet)
library(ape)
library(dplyr)
library(proxy)
library(stringr)
library(ggplot2)


vcf <- read.vcfR("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/trimmed_10bp_repeatmasked_vcf/lifted_vcf/lifted_all.vcf.gz")

# Extract DP (Depth) values from the VCF file
dp_matrix <- extract.gt(vcf, element = "DP", as.numeric = TRUE)

# Calculate average depth per variant (row-wise mean)
avg_depth_per_sample <- colMeans(dp_matrix, na.rm = TRUE)

write.csv(avg_depth_df, "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/trimmed_10bp_repeatmasked_vcf/lifted_all_depth.csv",
          row.names = FALSE)