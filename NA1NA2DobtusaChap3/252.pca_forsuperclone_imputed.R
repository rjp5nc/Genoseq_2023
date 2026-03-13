#ijob -A berglandlab -c2 -p standard --mem=80G
#module load gcc/11.4.0  openmpi/4.1.4 icu R/4.3.1
#R




  library(SNPRelate)
  library(SeqArray)
  library(gdsfmt)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(irlba)

# -------------- PARAMETERS ----------------
genofile.fn      <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/trimmed10bp_masked_usobtusa.gds"
meta.fn          <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/2022_2024seqmetadata20250811.csv"
samplestats.fn   <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/usdobtusa_samplestats.csv"
genomic_type.fn  <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/genomic_types.csv"

include_locations <- c("P66","P63","P58","P62","Gilmer")
MIN_DEPTH   <- 3
SNP_MAF     <- 0.05
SNP_MISS    <- 0.10
NUM_THREADS <- 10
SAVE_PLOTS  <- TRUE
OUTDIR      <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/"
SEED        <- 20251030
# ------------------------------------------


#Run on chr 3,7,8

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
set.seed(SEED)

# --- Load metadata ---
metadata     <- read.csv(meta.fn, header = TRUE, stringsAsFactors = FALSE)
samplestats  <- read.csv(samplestats.fn, header = TRUE, stringsAsFactors = FALSE)
genomic_type <- read.csv(genomic_type.fn, header = TRUE, stringsAsFactors = FALSE)

# --- Open GDS ---
gds <- seqOpen(genofile.fn)

# --- Select samples: location + depth ---
samples_to_keep <- metadata %>%
  inner_join(samplestats, by = c("Well" = "sampleId")) %>%
  filter(accuratelocation %in% include_locations, meanDepth > MIN_DEPTH) %>%
  pull(Well)
cat("# of selected samples:", length(samples_to_keep), "\n")

# --- SNP QC ---
snp.id <- snpgdsSelectSNP(gds,
                          sample.id = samples_to_keep,
                          autosome.only = FALSE,
                          remove.monosnp = TRUE,
                          maf = SNP_MAF,
                          missing.rate = SNP_MISS,
                          verbose = TRUE)
cat("# of selected SNPs:", length(snp.id), "\n")


# Retrieve chromosome info for all SNPs
all_chr <- seqGetData(gds, "chromosome")
all_snp_ids <- seqGetData(gds, "variant.id")

# Combine and subset only desired chromosomes (3, 7, 8)
snp_info <- data.frame(variant.id = all_snp_ids, chr = all_chr, stringsAsFactors = FALSE)
snp_info_filtered <- snp_info %>% filter(chr %in% c("JAACYE010000003.1", "JAACYE010000007.1", "JAACYE010000008.1"))

# Keep intersection of selected SNPs and desired chromosomes
snp.id <- intersect(snp.id, snp_info_filtered$variant.id)

cat("# of selected SNPs after QC and chromosome filter:", length(snp.id), "\n")



# --- Identify groups ---
groups_tbl <- metadata %>%
  inner_join(genomic_type, by = c("Well" = "CloneA")) %>%
  filter(Well %in% samples_to_keep) %>%
  select(Well, Group)

groups_tbl <- subset(groups_tbl, Group == "A")


# --- Helper: mean-impute + PCA ---
run_imputed_pca <- function(gds, sample_ids, snp_ids, n_pc = 10) {
  G <- snpgdsGetGeno(gds, sample.id = sample_ids, snp.id = snp_ids,
                     with.id = FALSE, snpfirstdim = TRUE)
  snp_means <- rowMeans(G, na.rm = TRUE)
  snp_means[!is.finite(snp_means)] <- 0
  snp_sds <- apply(G, 1, sd, na.rm = TRUE)
  snp_sds[!is.finite(snp_sds)] <- 0
  nas <- which(is.na(G), arr.ind = TRUE)
  if (nrow(nas) > 0) G[nas] <- snp_means[nas[,1]]
  keep <- snp_sds > 0
  Gz <- sweep(G[keep,,drop=FALSE], 1, snp_means[keep], "-")
  Gz <- sweep(Gz, 1, snp_sds[keep], "/")
  k <- min(n_pc, ncol(Gz)-1, nrow(Gz)-1)
  pca <- irlba::prcomp_irlba(t(Gz), n = k, center = FALSE, scale. = FALSE)
  list(sample.id = sample_ids, eigenvect = pca$x, eigenval = (pca$sdev)^2,
       kept_snps = sum(keep))
}

# --- Run PCA per group ---
var_summary <- data.frame(
  Group = character(),
  eigen_PC1 = numeric(), eigen_PC2 = numeric(), eigen_PC3 = numeric(),
  total_variance = numeric(),
  n_samples = integer(),
  kept_snps = integer(),
  stringsAsFactors = FALSE
)

pca_per_group_samples <- data.frame()
group_levels <- sort(unique(groups_tbl$Group))

for (g in group_levels) {
  cat("Running PCA (imputed) for Group =", g, "\n")
  samples.g <- groups_tbl %>% filter(Group == g) %>% pull(Well)
  n_samps <- length(samples.g)
  if (n_samps < 2) {
    cat("  Skipping (n < 2)\n")
    next
  }

  pca_g <- run_imputed_pca(gds, sample_ids = samples.g, snp_ids = snp.id, n_pc = 10)
  eigvals <- pca_g$eigenval
  eigen_PC1 <- eigvals[1]
  eigen_PC2 <- eigvals[2]
  eigen_PC3 <- eigvals[3]
  total_var <- sum(eigvals, na.rm = TRUE)

  var_summary <- rbind(var_summary,
    data.frame(Group = g,
               eigen_PC1 = eigen_PC1,
               eigen_PC2 = eigen_PC2,
               eigen_PC3 = eigen_PC3,
               total_variance = total_var,
               n_samples = n_samps,
               kept_snps = pca_g$kept_snps,
               stringsAsFactors = FALSE))

  pcs <- pca_g$eigenvect
  PC1 <- pcs[,1]
  PC2 <- if (ncol(pcs) >= 2) pcs[,2] else rep(0, length(PC1))
  PC3 <- if (ncol(pcs) >= 3) pcs[,3] else rep(0, length(PC1))

  temp <- data.frame(sample = pca_g$sample.id, PC1 = PC1, PC2 = PC2, PC3 = PC3,
                     Group = g, stringsAsFactors = FALSE)
  pca_per_group_samples <- rbind(pca_per_group_samples, temp)

  if (SAVE_PLOTS) {
    pfn <- file.path(OUTDIR, paste0("groupPCA_", g, "_PC1_PC2.png"))
    png(pfn, width = 2000, height = 1600, res = 200)
    plot(temp$PC1, temp$PC2, pch = 19,
         main = paste0("Group ", g, " PCA (imputed; n=", n_samps,
                       ", SNPs=", pca_g$kept_snps, ")"),
         xlab = paste0("PC1 (λ=", signif(eigen_PC1,3), ")", "PC1 (", round(100 * eigen_PC1 / total_var, 1), "% variance)"),
         ylab = paste0("PC2 (λ=", signif(eigen_PC2,3), ")", "PC2 (", round(100 * eigen_PC2 / total_var, 1), "% variance)"))
 #   text(temp$PC1, temp$PC2, labels = temp$sample, cex = 0.6, pos = 3)
    dev.off()
  }
}




# --- Save outputs ---
write.csv(var_summary,
          file.path(OUTDIR, "Group_PCA_variance_summary.csv"),
          row.names = FALSE)
write.csv(pca_per_group_samples,
          file.path(OUTDIR, "Group_PCA_scores_PC1-3.csv"),
          row.names = FALSE)

cat("Wrote group PCA results to:\n  - Group_PCA_variance_summary.csv\n  - Group_PCA_scores_PC1-3.csv\n")
seqClose(gds)
cat("Done.\n")
