# ijob -A berglandlab -c10 -p standard --mem=40G

#module load gcc/11.4.0 openmpi/4.1.4  R/4.3.1;R
### libraries
  library(data.table)
  library(foreach)
  library(SeqArray)
  library(gdsfmt)
library(SNPRelate)
library(ggplot2)
#  library(doMC)
 # registerDoMC(20)


genofile.fn <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/trimmed10bp_masked_usobtusa.gds"
gds <- seqOpen(genofile.fn)

metadata <- read.csv("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/2022_2024seqmetadata20250811.csv", header = TRUE)
samplestats <- read.csv("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/usdobtusa_samplestats.csv")
mitotypes <- read.csv("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/mito_types.csv")
genomic_types <- read.csv("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/genomic_types.csv")


unique(metadata$accuratelocation)

include_locations <- c("P66", "P63", "P58","P62", "Gilmer")  # replace with your names

# Keep only samples in these locations
samples_to_keep <- metadata$Well[metadata$accuratelocation %in% include_locations]
samples_to_keep <- samples_to_keep[samples_to_keep %in% genomic_types$CloneA]

seqSetFilter(gds, sample.id = samples_to_keep)

pca_result <- snpgdsPCA(gds, sample.id = samples_to_keep, num.thread=10, autosome.only=FALSE)

# Extract the PCA results
eigenvalues <- pca_result$eigenvect
eigenvectors <- pca_result$eigenvect

# Print the eigenvalues (variance explained by each principal component)
print(eigenvalues)

# Print the eigenvectors (coordinates for each individual)
print(head(eigenvectors))

pca_data <- data.frame(
  sample.id = pca_result$sample.id,  
  PC1 = pca_result$eigenvect[, 1],
  PC2 = pca_result$eigenvect[, 2],
  PC3 = pca_result$eigenvect[, 3],
  PC4 = pca_result$eigenvect[, 4],
  PC5 = pca_result$eigenvect[, 5],
  varprop = pca_result$varprop
)



write.csv(pca_data,
          "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/usobtusa_pca.csv",
          row.names = FALSE, quote = FALSE)


pca_merged <- merge(pca_data, metadata, by.x = "sample.id", by.y = "Well", all.x = TRUE)
pca_merged2 <- merge(pca_merged, samplestats, by.x = "sample.id", by.y = "sampleId", all.x = TRUE)



p <- ggplot(pca_merged2, aes(x = PC1, y = PC2, label = sample.id, col=accuratelocation)) +
  geom_point() +
  theme_bw() +
  labs(
    title = "PCA of USobtusa Samples",
    x = paste0("PC1 14.2%"),
    y = paste0("PC2 8.1%")
  )

# Save as PNG
ggsave("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/usobtusa_pca_plot.png",
       plot = p, width = 7, height = 6, dpi = 300)




p <- ggplot(pca_merged2, aes(x = PC2, y = PC3, label = sample.id, col=accuratelocation)) +
  geom_point() +
  theme_bw() +
  labs(
    title = "PCA of USobtusa Samples",
    x = paste0("PC2 8.1%"),
    y = paste0("PC3 6.9%")
  )

ggsave("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/usobtusa_pca2_pca3_plot.png",
       plot = p, width = 7, height = 6, dpi = 300)


p <- ggplot(pca_merged2, aes(x = PC1, y = PC2, label = sample.id, col=missingRate)) +
  geom_point() +
  theme_bw() +
  labs(
    title = "PCA of USobtusa Samples",
    x = paste0("PC1 14.2%"),
    y = paste0("PC2 8.1%")
  )

# Save as PNG
ggsave("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/usobtusa_pca_plot_missingrate.png",
       plot = p, width = 7, height = 6, dpi = 300)








p <- ggplot(pca_merged2, aes(x = PC1, y = PC2, label = sample.id, col=meanDepth)) +
  geom_point() +
  theme_bw() +
  labs(
    title = "PCA of USobtusa Samples",
    x = paste0("PC1 14.2%"),
    y = paste0("PC2 8.1%")
  )

# Save as PNG
ggsave("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/usobtusa_pca_plot_depth.png",
       plot = p, width = 7, height = 6, dpi = 300)





p <- ggplot(pca_merged2, aes(x = PC1, y = PC2, label = sample.id, shape=accuratelocation,col=date)) +
  geom_point() +
  theme_bw() +
  labs(
    title = "PCA of USobtusa Samples",
    x = paste0("PC1 14.2%"),
    y = paste0("PC2 8.1%")
  )

# Save as PNG
ggsave("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/usobtusa_pca_plot_date.png",
       plot = p, width = 7, height = 6, dpi = 300)




pca_merged3 <- merge(pca_merged2, mitotypes, by.x = "sample.id", by.y = "CloneA", all.x = TRUE)


p <- ggplot(pca_merged3, aes(x = PC1, y = PC2, label = sample.id,col=Group)) +
  geom_point() +
  theme_bw() +
  labs(
    title = "PCA of USobtusa Samples",
    x = paste0("PC1 14.2%"),
    y = paste0("PC2 8.1%")
  )

# Save as PNG
ggsave("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/mitotype_pca.png",
       plot = p, width = 7, height = 6, dpi = 300)





p <- ggplot(pca_merged3, aes(x = meanDepth, y = missingRate,col=accuratelocation)) +
  geom_point() +
  theme_bw() +
  labs(
    title = "Depth by Missing Rate per sample",
    x = paste0("Depth"),
    y = paste0("Missingrate")
  )

# Save as PNG
ggsave("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/missingratexdepth.png",
       plot = p, width = 7, height = 6, dpi = 300)