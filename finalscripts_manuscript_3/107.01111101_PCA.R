#ijob -A berglandlab -c2 -p standard --mem=40G
#module load gcc/11.4.0  openmpi/4.1.4 icu R/4.3.1
#R

library(tidyverse)
library(SNPRelate)
library(gdsfmt)
library(ggplot2)
library(data.table)

suppressPackageStartupMessages(library(ggrepel))

# ---- inputs ----
vcf_in    <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/trimmed10bp_filtered_hetregions_Gilmer_12.subset_regions.vcf.gz"

out_prefix <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/Gilmer_signature_01111101_fromHetRegions"

# ---- 1) get samples with signature 01111101 ----
keep_samples <- fread("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/Gilmer_signature_01111101_samples.txt")[[1]]

stopifnot(length(keep_samples) > 1)  # need at least 2 for PCA

list_file <- paste0(out_prefix, "_samples.txt")
writeLines(keep_samples, list_file)
cat("Saved sample list:", list_file, "\n")

# ---- 2) VCF -> SNP GDS ----
# gds_out <- paste0(out_prefix, ".gds")
# if (file.exists(gds_out)) file.remove(gds_out)

# snpgdsVCF2GDS(vcf_in, gds_out, method = "biallelic.only", snpfirstdim = TRUE)



# ---- 3) open SNP GDS and subset to your samples ----
g <- snpgdsOpen(gds_out)
all_samples <- read.gdsn(index.gdsn(g, "sample.id"))
keep <- intersect(all_samples, keep_samples)
if (length(keep) < 2) {
  snpgdsClose(g)
  stop(paste("After intersect, found", length(keep), "samples. Need at least 2."))
}
cat("Keeping", length(keep), "samples\n")

# optional but recommended: drop monomorphic sites
snp_set <- snpgdsSelectSNP(
  g,
  sample.id        = keep,
  missing.rate     = 1,     # keep all loci, let PCA impute missing
  maf              = 0,     # do not drop by MAF
  remove.monosnp   = TRUE,   # drop monomorphic among kept samples
  autosome.only = FALSE
)

cat("Selected", length(snp_set), "SNPs after removing monomorphic sites\n")

# ---- 4) PCA (SNPRelate mean-imputes missing) ----
pca <- snpgdsPCA(
  g,
  sample.id     = keep,
  snp.id        = snp_set,
  autosome.only = FALSE,  # contigs not standard autosomes
  missing.rate  = NaN,    # keep all selected sites
  num.thread    = 8,
  eigen.cnt     = 10
)
snpgdsClose(g)

# ---- 5) save outputs ----
pc_df <- tibble(
  sample = pca$sample.id,
  PC1 = pca$eigenvect[, 1],
  PC2 = pca$eigenvect[, 2],
  PC3 = pca$eigenvect[, 3]
)

out_csv <- paste0(out_prefix, "_PCA_scores.csv")
write.csv(pc_df, out_csv, row.names = FALSE)

# quick plot
pct <- round(100 * pca$varprop[1:2], 2)
p <- ggplot(pc_df, aes(PC1, PC2, label = sample)) +
  geom_point(size = 3) +
  ggrepel::geom_text_repel(size = 3, max.overlaps = 50) +
  theme_bw() +
  labs(
    title = "PCA on het-region SNPs (signature 01111101)",
    x = paste0("PC1 (", pct[1], "%)"),
    y = paste0("PC2 (", pct[2], "%)")
  )

out_png <- paste0(out_prefix, "_PCA.png")
ggsave(out_png, p, width = 7, height = 5, dpi = 300)

cat("Wrote:", out_csv, "\n")
cat("Wrote:", out_png, "\n")


