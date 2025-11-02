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

# Load metadata
metadata <- read.csv("/project/berglandlab/Robert/UKSequencing2022_2024/old_stuff/2022_2024seqmetadata20250131.csv", header = TRUE)

gilmer_good_samples <- subset(metadata, 
                              accuratelocation != "Gilmer" & !(clone %in% c("BLANK", "Blank")))$Well

# Get all sample IDs from the GDS
all_samples <- seqGetData(genofile, "sample.id")

# Intersect with GDS sample IDs to avoid mismatches
final_samples <- intersect(all_samples, gilmer_good_samples)

# Apply sample filter to GDS
seqSetFilter(genofile, sample.id = final_samples, verbose = TRUE)


# ===== Missing Rate per Sample =====
miss_sample <- seqMissing(genofile, per.variant=FALSE)


png(file.path(output_dir, "missing_rate_per_sample_gilmer.png"), width=800, height=600)
hist(miss_sample, breaks=40, main="Missing Rate per Sample", xlab="Missing Rate")
abline(v=0.1, col="red", lty=2)
dev.off()

# ===== Missing Rate per Variant =====
miss_variant <- seqMissing(genofile, per.variant=TRUE)

png(file.path(output_dir, "missing_rate_per_variant_gilmer.png"), width=800, height=600)
hist(miss_variant, breaks=40, main="Missing Rate per Variant", xlab="Missing Rate")
abline(v=0.1, col="red", lty=2)
dev.off()

# ===== Filter by MAF & Missingness =====
maf <- seqAlleleFreq(genofile)
maf_filter <- (maf >= 0.05 & maf <= 0.95)
miss_filter <- (miss_variant <= 0.10)
var_filter <- which(maf_filter & miss_filter)
seqSetFilter(genofile, sample.id = final_samples, variant.id = var_filter)

# ===== PCA =====
pca <- snpgdsPCA(genofile, sample.id = final_samples, autosome.only = FALSE)
pc.percent <- pca$varprop * 100

# Save PCA plot
png(file.path(output_dir, "pca_plot_not_gilmer.png"), width=800, height=600)
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
seqSetFilter(genofile, sample.id = final_samples, variant.id = subset_snp_ids, verbose = FALSE)

sample_ids <- seqGetData(genofile, "sample.id")  # Extract sample IDs

# Compute IBS matrix on a SNP subset (you already did this)
ibd_dist <- snpgdsIBS(genofile, autosome.only = FALSE, verbose = TRUE)
ibs_matrix <- 1 - ibd_dist$ibs  # Convert similarity to distance

# Build tree
tree <- nj(as.dist(ibs_matrix))

# Assign sample names to tip labels
tree$tip.label <- sample_ids

# Plot and save tree
png(file.path(output_dir, "nj_tree_not_gilmer.png"), width=800, height=600)
plot(tree, main = "Neighbor-Joining Tree with Sample Names")
dev.off()

# Save tree object and sample names
save(tree, file = file.path(output_dir, "nj_tree_not_gilmer.RData"))
write.csv(data.frame(Sample = tree$tip.label),
          file.path(output_dir, "nj_tree_sample_names_not_gilmer.csv"),
          row.names = FALSE)



# ===== Export Data =====
write.csv(data.frame(Sample = pca$sample.id,
                     PC1 = pca$eigenvect[,1],
                     PC2 = pca$eigenvect[,2]),
          file.path(output_dir, "seqarray_pca_not_gilmer.csv"), row.names = FALSE)



ibs <- snpgdsIBS(genofile, sample.id = final_samples, autosome.only = FALSE)
ibs_matrix <- ibs$ibs  # this is the numeric similarity matrix
rownames(ibs_matrix) <- colnames(ibs_matrix) <- ibs$sample.id

write.csv(ibs_matrix, file.path(output_dir, "ibs_matrix_not_gilmer.csv"))






set.seed(100)

ibd <- snpgdsIBDMoM(genofile, sample.id=gilmer_good_samples, snp.id=snp.id,
    maf=0.05, missing.rate=0.05, num.thread=2, verbose = TRUE, autosome.only = FALSE)

relatedness <- snpgdsIBDSelection(ibd)
write.csv(relatedness, file.path(output_dir, "estimated_relatedness_not_gilmer.csv"), row.names = FALSE)

# Reset filter to full variant set (if needed later)
seqResetFilter(genofile, verbose = FALSE)
# Save IBD data

seqClose(genofile)

cat("✅ QC complete. Plots and results saved in:", output_dir, "\n")





