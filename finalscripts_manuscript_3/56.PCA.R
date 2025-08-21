# ijob -A berglandlab -c10 -p standard --mem=40G
### module load gcc/7.1.0  openmpi/3.1.4 R/4.1.1; R


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

pca_result <- snpgdsPCA(gds, num.thread=10, autosome.only=FALSE)

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


p <- ggplot(pca_merged, aes(x = PC1, y = PC2, label = sample.id, col=accuratelocation)) +
  geom_point() +
  theme_bw() +
  labs(
    title = "PCA of USobtusa Samples",
    x = paste0("PC1"),
    y = paste0("PC2")
  )

# Save as PNG
ggsave("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/usobtusa_pca_plot.png",
       plot = p, width = 7, height = 6, dpi = 300)