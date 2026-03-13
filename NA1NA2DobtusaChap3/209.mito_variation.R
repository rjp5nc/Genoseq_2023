#ijob -A berglandlab -c10 -p standard --mem=40G
#module load gcc/11.4.0  openmpi/4.1.4 icu R/4.3.1; R

#14642bp

suppressPackageStartupMessages({
  library(SeqArray)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(ape)
  library(ggtree)
  library(ggplot2)
  library(scales)
})

# ----------------------------
# Paths
# ----------------------------
gds_path <- "/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/usdobtusa_mito_allsites_all.haploid.annotated2.gds"
mitotype_csv <- "/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/mito_types.csv"
out_dir <- "/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/unique_snps_by_group"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# Parameters
# ----------------------------
dp_cutoff      <- 20
miss_cutoff    <- 0.05
maf_cutoff     <- 0.00
impute_missing <- TRUE
root_outgroup  <- "SRR5012393"  # set NULL to skip

# ----------------------------
# Read mitotypes (robust)
# Expects 2 columns: sample, mitotype/group
# ----------------------------
mito_raw <- read_csv(mitotype_csv, show_col_types = FALSE)

stopifnot(ncol(mito_raw) >= 2)

mito <- mito_raw %>%
  as.data.frame() %>%                 # in case it's a matrix
  tibble::as_tibble() %>%
  { setNames(., c("sample", "mito_type", names(.)[-(1:2)])) } %>%  # rename first 2
  transmute(
    sample    = str_trim(as.character(sample)),
    mito_type = str_trim(as.character(mito_type))
  ) %>%
  filter(nzchar(sample)) %>%
  mutate(mito_plot = ifelse(is.na(mito_type) | mito_type == "", "Unknown", mito_type)) %>%
  distinct(sample, .keep_all = TRUE)


# ----------------------------
# Open GDS (only once)
# ----------------------------
options(mc.cores = 1)
Sys.setenv(OMP_NUM_THREADS = "1")

g <- seqOpen(gds_path)
on.exit(seqClose(g), add = TRUE)
seqResetFilter(g)

# ----------------------------
# Sample filter by mean DP
# ----------------------------
DP <- seqGetData(g, "annotation/format/DP")
if (is.list(DP)) DP <- do.call(cbind, DP)

samps_all <- seqGetData(g, "sample.id")

mean_dp <- if (ncol(DP) == length(samps_all)) colMeans(DP, na.rm = TRUE) else rowMeans(DP, na.rm = TRUE)
keep_samples <- samps_all[!is.na(mean_dp) & mean_dp > dp_cutoff]

seqResetFilter(g)
seqSetFilter(g, sample.id = keep_samples, action = "set")
samps <- seqGetData(g, "sample.id")
message("Samples kept after DP filter: ", length(samps))

# ----------------------------
# Genotypes (haploid) -> matrix variants x samples
# ----------------------------
GT_raw <- seqGetData(g, "genotype")
variant_ids <- seqGetData(g, "variant.id")

if (length(dim(GT_raw)) == 3L) {
  d <- dim(GT_raw)
  # ploidy x variant x sample
  if (d[2] == length(variant_ids) && d[3] == length(samps)) {
    GT <- matrix(GT_raw[1, , ], nrow = d[2], ncol = d[3])
  # ploidy x sample x variant
  } else if (d[2] == length(samps) && d[3] == length(variant_ids)) {
    GT <- matrix(GT_raw[1, , ], nrow = d[2], ncol = d[3])
    GT <- t(GT)
  } else {
    stop("Unexpected genotype array shape: ", paste(d, collapse="x"))
  }
} else if (is.matrix(GT_raw)) {
  GT <- GT_raw
  if (ncol(GT) != length(samps) && nrow(GT) == length(samps)) GT <- t(GT)
} else stop("Unexpected genotype object type.")

stopifnot(ncol(GT) == length(samps))

# clean coding
GT[is.na(GT) | GT >= 3L] <- NA_integer_
GT[!is.na(GT) & GT > 1L] <- 1L

# ----------------------------
# Variant filters: polymorphic + missingness (+optional MAF)
# ----------------------------
is_poly <- apply(GT, 1, function(z) {
  u <- unique(z[!is.na(z)])
  length(u) > 1
})
GT <- GT[is_poly, , drop = FALSE]

miss_rate <- rowMeans(is.na(GT))
GT <- GT[miss_rate <= miss_cutoff, , drop = FALSE]

if (maf_cutoff > 0) {
  af  <- rowMeans(GT, na.rm = TRUE)
  maf <- pmin(af, 1 - af)
  GT  <- GT[maf >= maf_cutoff, , drop = FALSE]
}

stopifnot(nrow(GT) >= 2, ncol(GT) == length(samps))
message("Variants kept for PCA/tree: ", nrow(GT))

# ----------------------------
# PCA (impute missing only here)
# ----------------------------
GT_for_pca <- GT
if (impute_missing && anyNA(GT_for_pca)) {
  mu <- rowMeans(GT_for_pca, na.rm = TRUE)
  na_idx <- which(is.na(GT_for_pca), arr.ind = TRUE)
  GT_for_pca[na_idx] <- mu[na_idx[, 1]]
}

X <- t(GT_for_pca)  # samples x variants

pca <- prcomp(X, center = TRUE, scale. = FALSE)
varprop <- (pca$sdev^2) / sum(pca$sdev^2)

pca_df <- as.data.frame(pca$x[, 1:min(5, ncol(pca$x)), drop = FALSE]) %>%
  rownames_to_column("sample") %>%
  left_join(mito %>% select(sample, mito_plot), by = "sample") %>%
  mutate(mito_plot = ifelse(is.na(mito_plot), "Unknown", mito_plot))

write_csv(pca_df, file.path(out_dir, "pca_mtDNA_filtered.csv"))

p_pca <- ggplot(pca_df, aes(PC1, PC2, color = mito_plot)) +
  geom_point(size = 2) +
  theme_bw() +
  labs(
    title = "PCA (mtDNA, polymorphic sites, missingness-filtered)",
    x = paste0("PC1 (", round(varprop[1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(varprop[2] * 100, 1), "%)"),
    color = "Mitotype"
  )
ggsave(file.path(out_dir, "PCA_PC1_PC2.png"), p_pca, width = 7, height = 6, dpi = 300)

# ----------------------------
# IBS distance + NJ tree (NO IMPUTATION)
# ----------------------------
n <- ncol(GT)
IBS <- matrix(NA_real_, n, n, dimnames = list(samps, samps))
for (i in seq_len(n)) {
  IBS[i, i] <- 1
  if (i < n) {
    for (j in (i + 1L):n) {
      ok <- !(is.na(GT[, i]) | is.na(GT[, j]))
      denom <- sum(ok)
      v <- if (denom == 0) NA_real_ else sum(GT[ok, i] == GT[ok, j]) / denom
      IBS[i, j] <- IBS[j, i] <- v
    }
  }
}
D <- 1 - IBS
if (anyNA(D)) D[is.na(D)] <- max(D, na.rm = TRUE)

tree <- nj(as.dist(D))
tip_meta <- tibble(label = tree$tip.label) %>%
  left_join(mito %>% select(sample, mito_plot), by = c("label" = "sample")) %>%
  mutate(mito_plot = ifelse(is.na(mito_plot), "Unknown", mito_plot))

mt_levels <- sort(unique(tip_meta$mito_plot))
pal <- setNames(hue_pal()(length(mt_levels)), mt_levels)

p_rect <- ggtree(tree) %<+% tip_meta +
  geom_tiplab(aes(color = mito_plot), size = 2) +
  scale_color_manual(values = pal, name = "Mitotype") +
  theme_tree2() +
  labs(title = "NJ tree (IBS distance, mtDNA polymorphic sites)")
ggsave(file.path(out_dir, "NJ_IBS_rect.png"), p_rect, width = 10, height = 10, dpi = 300)

p_circ <- ggtree(tree, layout = "circular") %<+% tip_meta +
  geom_tiplab(aes(color = mito_plot), size = 2) +
  scale_color_manual(values = pal, name = "Mitotype") +
  theme_tree2() +
  labs(title = "NJ tree (IBS distance, mtDNA polymorphic sites)")
ggsave(file.path(out_dir, "NJ_IBS_circular.png"), p_circ, width = 10, height = 10, dpi = 300)
# ----------------------------
# Optional rooting
# ----------------------------

root_outgroup <- "RobertUK_B10"

if (!is.null(root_outgroup) && root_outgroup %in% tree$tip.label) {
  tree_rooted <- root(tree, outgroup = root_outgroup, resolve.root = TRUE)

  p_rect_r <- ggtree(tree_rooted) %<+% tip_meta +
    geom_tiplab(aes(color = mito_plot), size = 2) +
    scale_color_manual(values = pal, name = "Mitotype") +
    theme_tree2() +
    labs(title = paste0("NJ tree rooted by ", root_outgroup))
  ggsave(file.path(out_dir, "NJ_IBS_rect_rooted.png"), p_rect_r, width = 10, height = 10, dpi = 300)

  p_circ_r <- ggtree(tree_rooted, layout = "circular") %<+% tip_meta +
    geom_tiplab(aes(color = mito_plot), size = 2) +
    scale_color_manual(values = pal, name = "Mitotype") +
    theme_tree2() +
    labs(title = paste0("NJ tree rooted by ", root_outgroup))
  ggsave(file.path(out_dir, "NJ_IBS_circular_rooted.png"), p_circ_r, width = 10, height = 10, dpi = 300)
} else {
  message("Rooting skipped (outgroup not found or root_outgroup is NULL).")
}