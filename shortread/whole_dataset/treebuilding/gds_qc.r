# Load libraries
library(SeqArray)
library(SNPRelate)
library(gdsfmt)
library(ggplot2)
library(ape)

# Set your file path
gds_file <- "/scratch/rjp5nc/UK2022_2024/gds_USobtusa_temp/raw_usobtusa.gds.gds.gz"

output_dir <- "/scratch/rjp5nc/UK2022_2024/gds_USobtusa_temp/"

# Open GDS file
genofile <- seqOpen(gds_file)

# ===== Basic Summary =====
seqSummary(genofile)

# ===== Missing Rate per Sample =====
miss_sample <- seqMissing(genofile, per.variant=FALSE)


png(file.path(output_dir, "missing_rate_per_sample.png"), width=800, height=600)
hist(miss_sample, breaks=40, main="Missing Rate per Sample", xlab="Missing Rate")
abline(v=0.1, col="red", lty=2)
dev.off()

# ===== Missing Rate per Variant =====
miss_variant <- seqMissing(genofile, per.variant=TRUE)

png(file.path(output_dir, "missing_rate_per_variant.png"), width=800, height=600)
hist(miss_variant, breaks=40, main="Missing Rate per Variant", xlab="Missing Rate")
abline(v=0.1, col="red", lty=2)
dev.off()

# ===== Filter by MAF & Missingness =====
maf <- seqAlleleFreq(genofile)
maf_filter <- (maf >= 0.05 & maf <= 0.95)
miss_filter <- (miss_variant <= 0.10)
var_filter <- which(maf_filter & miss_filter)
seqSetFilter(genofile, variant.id = var_filter)

# ===== PCA =====
pca <- snpgdsPCA(genofile, autosome.only = FALSE)
pc.percent <- pca$varprop * 100

# Save PCA plot
png(file.path(output_dir, "pca_plot.png"), width=800, height=600)
plot(pca$eigenvect[,1], pca$eigenvect[,2],
     col = "blue", pch = 19,
     xlab = paste0("PC1 (", round(pc.percent[1], 2), "%)"),
     ylab = paste0("PC2 (", round(pc.percent[2], 2), "%)"),
     main = "PCA of Sequenced Samples")
dev.off()

# ===== Random SNP Subset for IBS and Tree =====

# Select a random subset of 10,000 SNPs that passed previous filters
set.seed(123)  # for reproducibility
snp_ids <- seqGetData(genofile, "variant.id")
subset_snp_ids <- sample(snp_ids, size = min(10000, length(snp_ids)))

# Temporarily apply filter
seqSetFilter(genofile, variant.id = subset_snp_ids, verbose = FALSE)

# Compute IBS and build tree
ibd_dist <- snpgdsIBS(genofile, autosome.only = FALSE, verbose = TRUE)
ibs_matrix <- 1 - ibd_dist$ibs  # Convert similarity to distance
tree <- nj(as.dist(ibs_matrix))

# Plot and save NJ tree
png(file.path(output_dir, "nj_tree.png"), width=800, height=600)
plot(tree, main = "Neighbor-Joining Tree from IBS Distance (10K SNPs)")
dev.off()

# Reset filter to full variant set (if needed later)
seqResetFilter(genofile, verbose = FALSE)


# ===== Export Data =====
write.csv(data.frame(Sample = pca$sample.id,
                     PC1 = pca$eigenvect[,1],
                     PC2 = pca$eigenvect[,2]),
          file.path(output_dir, "seqarray_pca.csv"), row.names = FALSE)

# Optionally save tree object
save(tree, file = file.path(output_dir, "nj_tree.RData"))


seqClose(genofile)

cat("✅ QC complete. Plots and results saved in:", output_dir, "\n")
