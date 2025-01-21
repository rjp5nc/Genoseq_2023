# ijob -A berglandlab -c20 -p standard --mem=40G
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



genofile.fn <- "/project/berglandlab/Robert/UKSequencing2022_2024/usftp21.novogene.com/01.RawData/Bams/vcf/2022seq.concat.Removereps.gds.gz"
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
                       PC2 = pca_result$eigenvectors[, 2])

# Plot PCA (PC1 vs PC2)

write.csv(pca_data,"/project/berglandlab/Robert/2022seqPCA.csv")
