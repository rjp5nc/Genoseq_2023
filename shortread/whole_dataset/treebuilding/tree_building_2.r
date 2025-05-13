library(vcfR)
library(adegenet)
library(ape)
library(proxy)
library(ggplot2)
library(stringr)

# Read VCF file
vcf <- read.vcfR("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/raw_vcf/raw_usobtusa_vcf.vcf.gz")

# Extract depth (DP) values
dp_matrix <- extract.gt(vcf, element = "DP", as.numeric = TRUE)

# Calculate average depth per variant (row-wise mean)
avg_depth_per_variant <- rowMeans(dp_matrix, na.rm = TRUE)

# Filter: Keep only variants with average depth > 4
variants_to_keep <- which(avg_depth_per_variant > 4)
vcf_filtered <- vcf[variants_to_keep, ]

# Extract genotype data from filtered VCF
genotypes <- extract.gt(vcf_filtered)

# Convert to genind object
genind_obj <- vcfR2genind(vcf_filtered)

# Load metadata
metadata <- read.csv("/project/berglandlab/Robert/UKSequencing2022_2024/old_stuff/2022_2024seqmetadata20250131.csv", header = TRUE)

# Clean metadata (remove blanks)
metadata_with_clone <- subset(metadata, !clone %in% c("BLANK", "Blank"))

# Subset for Rockpool samples
samples_to_keep <- metadata_with_clone$Well[metadata_with_clone$location == "Rockpool" | metadata_with_clone$location == "Gilmer"]

# Subset genind object to Rockpool samples
genind_obj2 <- genind_obj[indNames(genind_obj) %in% samples_to_keep, ]

# Convert genind to matrix and compute distance
genind_matrix <- tab(genind_obj2, NA.method = "mean")
dist_matrix <- dist(genind_matrix)

# Build phylogenetic tree
tree <- njs(dist_matrix)
save(tree, file = "/scratch/rjp5nc/UK2022_2024/daptree/phylogenetic_treewgs4x.RData")

# Compute Jaccard similarity matrix
sim_matrix <- as.matrix(simil(genind_matrix, method = "Jaccard"))

# Save similarity matrix
write.csv(sim_matrix, "/scratch/rjp5nc/UK2022_2024/daptree/similarity_matrix4x.csv")

# Optional: Heatmap
heatmap(sim_matrix)

# Depth plot (optional part at the end, update if needed)
# If you still want to generate average depth histogram across samples, you'll need to calculate colMeans:
avg_depth_per_sample <- colMeans(dp_matrix, na.rm = TRUE)

depth_data <- data.frame(Sample = names(avg_depth_per_sample), 
                         Depth = as.numeric(avg_depth_per_sample), 
                         stringsAsFactors = FALSE)
depth_data$Sample <- str_remove(depth_data$Sample, ".sort.dedup")

p <- ggplot(depth_data, aes(x = Depth)) +
  geom_histogram(binwidth = 1, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Histogram of Read Depth", x = "Average Read Depth", y = "Frequency") +
  theme_bw()

ggsave("/scratch/rjp5nc/UK2022_2024/daptree/avg_depth_histogram.png", plot = p, width = 6, height = 4, dpi = 300)
