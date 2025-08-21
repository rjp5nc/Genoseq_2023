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


pca_result <- snpgdsPCA(gds, num.thread=10, autosome.only=FALSE)

# Extract the PCA results
eigenvalues <- pca_result$eigenvect
eigenvectors <- pca_result$eigenvect

# Print the eigenvalues (variance explained by each principal component)
print(eigenvalues)

# Print the eigenvectors (coordinates for each individual)
print(head(eigenvectors))

pc_percent <- round(pca_result$varprop * 100, 2)

pca_data <- data.frame(
  sample.id = pca_result$sample.id,  
  PC1 = pca_result$eigenvect[, 1],
  PC2 = pca_result$eigenvect[, 2],
  PC3 = pca_result$eigenvect[, 3],
  PC4 = pca_result$eigenvect[, 4],
  PC5 = pca_result$eigenvect[, 5]
)

colnames(pca_data) <- c(
  "sample.id",
  paste0("PC1 (", pc_percent[1], "%)"),
  paste0("PC2 (", pc_percent[2], "%)"),
  paste0("PC3 (", pc_percent[3], "%)"),
  paste0("PC4 (", pc_percent[4], "%)"),
  paste0("PC5 (", pc_percent[5], "%)")
)

write.csv(pca_data,
          "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/usobtusa_pca.csv",
          row.names = FALSE, quote = FALSE)


p <- ggplot(pca_data, aes(x = PC1, y = PC2, label = sample.id)) +
  geom_point(color = "steelblue", size = 2, alpha = 0.7) +
  theme_minimal(base_size = 14) +
  labs(
    title = "PCA of Daphnia Samples",
    x = paste0("PC1 (", pc_percent[1], "%)"),
    y = paste0("PC2 (", pc_percent[2], "%)")
  )

# Save as PNG
ggsave("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/usobtusa_pca_plot.png",
       plot = p, width = 7, height = 6, dpi = 300)