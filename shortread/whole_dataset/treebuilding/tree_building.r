
library(vcfR)
library(adegenet)
library(ape)
library(dplyr)
library(proxy)
library(stringr)
library(ggplot2)


vcf <- read.vcfR("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/raw_vcf/raw_usobtusa_vcf.vcf.gz")

# Extract DP (Depth) values from the VCF file
dp_matrix <- extract.gt(vcf, element = "DP", as.numeric = TRUE)

# Calculate average depth per variant (row-wise mean)
#avg_depth_per_sample <- colMeans(dp_matrix, na.rm = TRUE)

# Filter variants where the average depth is above 4

# Step 2: Extract Genotype Data
genotypes <- extract.gt(vcf)

# Step 3: Convert Genotype Data to Genind Object
genind_obj <- vcfR2genind(vcf)


metadata <- read.csv("/project/berglandlab/Robert/UKSequencing2022_2024/old_stuff/2022_2024seqmetadata20250131.csv", header = TRUE)

metadata_with_clone2 <- subset(metadata, clone != "BLANK")
metadata_with_clone3 <- subset(metadata_with_clone2, clone != "Blank")

# Filter metadata to include only samples from "Rockpool"
samples_to_keep <- metadata_with_clone3$Well[metadata_with_clone3$location == "Rockpool"]

head(metadata)


## Step 6: Subset Genind Object to Keep Only Selected Samples
genind_obj2 <- genind_obj[indNames(genind_obj) %in% samples_to_keep, ]

# Step 4: Compute Distance Matrix
genind_obj <- as.matrix(genind_obj)

dist_matrix <- dist(genind_obj)


# Step 5: Build Phylogenetic Tree (Neighbor-Joining)
tree <- njs(dist_matrix)
save(tree, file = "/scratch/rjp5nc/UK2022_2024/daptree/phylogenetic_treewgs.RData")



#sim_matrix <- 1 / (1 + (dist_matrix))  # Convert to similarity (scaled 0-1)

#genind_obj3<- tab(genind_obj2)
sim_matrix <- as.matrix(simil(genind_obj, method = "Jaccard"))

heatmap(sim_matrix)  # Quick visualization

#genind_data <- tab(genind_obj)  # Convert genind object to genotype table
sim_matrix <- simil(genind_data, method = "Jaccard")

# Optionally, convert to a numeric matrix for further processing
sim_matrix <- as.matrix(sim_matrix)

write.csv(sim_matrix, "/scratch/rjp5nc/UK2022_2024/daptree/similarity_matrix.csv")  # Save as CSV











depth_data <- data.frame(Sample = names(avg_depth_per_sample), 
                         Depth = gsub("\\.sort\\.dedup$", "", avg_depth_per_sample), 
                         stringsAsFactors = FALSE)
depth_data$Sample <- str_remove(depth_data$Sample, ".sort.dedup")
depth_data$Depth <- as.numeric(depth_data$Depth)

p <- ggplot(depth_data, aes(x = Depth)) +
  geom_histogram(binwidth = 1, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Histogram of Read Depth", x = "Average Read Depth", y = "Frequency") +
  theme_bw()

# Save plot as PNG
ggsave("/scratch/rjp5nc/UK2022_2024/allshortreads/forspecies/avg_depth_histogram.png", plot = p, width = 6, height = 4, dpi = 300)
