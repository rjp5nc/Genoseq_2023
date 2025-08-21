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
eigenvalues <- pca_result$eigenvalues
eigenvectors <- pca_result$eigenvectors

# Print the eigenvalues (variance explained by each principal component)
print(eigenvalues)

# Print the eigenvectors (coordinates for each individual)
print(head(eigenvectors))

pca_data <- data.frame(PC1 = pca_result$eigenvectors[, 1], 
                       PC2 = pca_result$eigenvectors[, 2], 
                       PC3 = pca_result$eigenvectors[, 3], 
                       PC4 = pca_result$eigenvectors[, 4], 
                       PC5 = pca_result$eigenvectors[, 5])

# Plot PCA (PC1 vs PC2)

write.csv(pca_data,"/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/usobtusa_pca.csv")
