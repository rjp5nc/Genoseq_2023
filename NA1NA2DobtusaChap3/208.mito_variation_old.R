#ijob -A berglandlab -c10 -p standard --mem=40G
#module load gcc/11.4.0  openmpi/4.1.4 icu R/4.3.1
#R

#14642bp



library(SeqArray)
  library(SeqVarTools)
  library(gdsfmt)
  library(ggplot2)
  library(dplyr)
library(GenomicRanges)
library(rtracklayer)



gds_file <- "/scratch/rjp5nc/UK2022_2024/allsites_mito/usdobtusa_mito_allsites.haploid.gds"
out_dir  <- "/scratch/rjp5nc/UK2022_2024/allsites_mito/plots_diff_per_sample"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

gds <- seqOpen(gds_file)

# ---- Basic metadata ----
samples <- seqGetData(gds, "sample.id")
chrom   <- seqGetData(gds, "chromosome")
pos     <- seqGetData(gds, "position")
n_samp  <- length(samples)
n_var   <- length(pos)

message(sprintf("Loaded: %s | %d samples | %d sites", basename(gds_file), n_samp, n_var))

# ---- Helper: fetch per-sample depth matrix ----
get_depth_matrix <- function(gds) {
  has_dp <- gdsfmt::index.gdsn(gds, "annotation/format/DP", silent = TRUE)
  if (!is.null(has_dp)) {
    dp <- seqGetData(gds, "annotation/format/DP")
    # DP can come as a list or array; coerce to matrix [samples x variants] if needed
    if (is.list(dp) && !is.null(dp$data)) dp <- dp$data
    dp <- as.matrix(dp)
    # Try to standardize orientation to [samples x variants]
    if (nrow(dp) == length(pos) && ncol(dp) == length(samples)) {
      dp <- t(dp)
    }
    return(dp)
  }
  # Fallback to AD (allelic depths): sum across alleles to approximate DP
  has_ad <- gdsfmt::index.gdsn(gds, "annotation/format/AD", silent = TRUE)
  if (!is.null(has_ad)) {
    ad <- seqGetData(gds, "annotation/format/AD")
    # AD is often [allele x sample x variant]; sum over allele
    ad <- apply(ad, c(2,3), function(x) sum(x, na.rm = TRUE))
    return(ad) # [samples x variants]
  }
  stop("No per-sample depth found (neither FORMAT/DP nor FORMAT/AD present).")
}

# ---- Compute mean depth per sample, filter ≥ 20 ----
dp <- get_depth_matrix(gds)
stopifnot(nrow(dp) == n_samp, ncol(dp) == n_var)

mean_dp <- rowMeans(dp, na.rm = TRUE)
keep_idx <- which(mean_dp >= 20)
keep_samples <- samples[keep_idx]
message(sprintf("Samples passing mean DP ≥ 20: %d / %d", length(keep_samples), length(samples)))

if (length(keep_samples) == 0) {
  stop("No samples meet mean depth ≥ 20.")
}

# ---- Get haploid genotypes and build difference indicator ----
# Genotype array is usually [ploidy x sample x variant]; for haploid, ploidy = 1 (0=ref, 1=alt, 2+=other/missing codes)
gt <- seqGetData(gds, "genotype")        # array
dim_gt <- dim(gt)
if (length(dim_gt) != 3) stop("Unexpected genotype dimensions.")
if (dim_gt[1] != 1) warning("Ploidy != 1 detected; proceeding but assuming first allele represents haploid call.")

# Extract to matrix [sample x variant]
gt_mat <- gt[1, , , drop = FALSE]
gt_mat <- matrix(gt_mat, nrow = n_samp, ncol = n_var, byrow = FALSE)

# Difference from reference: 1 if non-ref (==1) else 0; treat negative or 3 as missing
is_miss <- is.na(gt_mat) | gt_mat < 0
diff_mat <- ifelse(!is_miss & gt_mat != 0, 1L, 0L)  # ref=0, alt=1 → 1 indicates difference

# ---- Plot per-sample genome-wide differences ----
# We facet by chromosome; to keep files reasonably sized, we build per-sample data.frame on the fly
chrom_factor <- factor(chrom, levels = unique(chrom))

for (i in keep_idx) {
  sname <- samples[i]
  diffs <- diff_mat[i, ]
  df <- data.frame(
    chromosome = chrom_factor,
    position   = pos,
    diff_from_ref = diffs
  )
  # Drop missing positions if any (should already be 0/1)
  df <- df[!is.na(df$diff_from_ref), ]

  p <- ggplot(df, aes(x = position, y = diff_from_ref)) +
    geom_point(size = 0.2, alpha = 0.4) +
    facet_wrap(~ chromosome, scales = "free_x", ncol = 1) +
    scale_y_continuous(breaks = c(0, 1), limits = c(-0.05, 1.05)) +
    labs(
      title = paste0("Genome-wide differences from reference (haploid) — ", sname),
      x = "Genomic position",
      y = "Difference from ref (0/1)"
    ) +
    theme_bw(base_size = 10) +
    theme(
      strip.background = element_rect(fill = "grey90"),
      panel.grid.minor = element_blank()
    )

  out_png <- file.path(out_dir, sprintf("diff_from_ref_%s.png", gsub("[^A-Za-z0-9._-]", "_", sname)))
  ggsave(out_png, p, width = 10, height = max(3.5, length(levels(chrom_factor)) * 1.8), dpi = 250)
  message(sprintf("Wrote: %s", out_png))
}

seqClose(gds)
message("Done.")
















suppressPackageStartupMessages({
  library(SeqArray)
  library(SeqVarTools)
  library(gdsfmt)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(purrr)
})

gds_file <- "/scratch/rjp5nc/UK2022_2024/allsites_mito/usdobtusa_mito_allsites.haploid.gds"
out_dir  <- "/scratch/rjp5nc/UK2022_2024/allsites_mito/"
out_pdf  <- file.path(out_dir, "pairwise_differences_all_pairs.pdf")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
gds <- seqOpen(gds_file)

# ---- metadata ----
samples <- seqGetData(gds, "sample.id")
chrom   <- seqGetData(gds, "chromosome")
pos     <- seqGetData(gds, "position")
n_samp  <- length(samples)
n_var   <- length(pos)
chrom_levels <- unique(chrom)

message(sprintf("Loaded %s | %d samples | %d sites",
                basename(gds_file), n_samp, n_var))

# ---- depth helper ----
get_depth_matrix <- function(gds, pos_len, samp_len, chrom=NULL, pos=NULL) {
  has_dp <- gdsfmt::index.gdsn(gds, "annotation/format/DP", silent = TRUE)
  if (!is.null(has_dp)) {
    dp <- seqGetData(gds, "annotation/format/DP")
    if (is.list(dp) && !is.null(dp$data)) dp <- dp$data
    dp <- as.matrix(dp)
    if (nrow(dp) == pos_len && ncol(dp) == samp_len) dp <- t(dp)
    return(dp)
  }
  has_ad <- gdsfmt::index.gdsn(gds, "annotation/format/AD", silent = TRUE)
  if (!is.null(has_ad)) {
    ad <- seqGetData(gds, "annotation/format/AD")  # [allele x sample x variant]
    dp <- apply(ad, c(2,3), function(x) sum(x, na.rm = TRUE))
    return(dp)
  }
  stop("No per-sample depth found (neither FORMAT/DP nor FORMAT/AD).")
}

# ---- filter samples by mean DP >= 20 ----
dp <- get_depth_matrix(gds, pos_len = n_var, samp_len = n_samp)
mean_dp <- rowMeans(dp, na.rm = TRUE)
keep_idx <- which(mean_dp >= 20)
keep_samples <- samples[keep_idx]
message(sprintf("Samples passing mean DP ≥ 20: %d / %d", length(keep_samples), n_samp))
if (!length(keep_samples)) stop("No samples meet mean depth ≥ 20.")

# ---- prep pairs ----
pair_mat <- t(combn(keep_samples, 2))
pair_labels <- paste(pair_mat[,1], "vs", pair_mat[,2])
n_pairs <- nrow(pair_mat)
message(sprintf("Total pairwise comparisons: %d", n_pairs))

# ---- iterate variants by chunk, collect mismatches sparsely ----
chunk_size <- 5000L
variant_ids <- seq_len(n_var)

collect_pairs <- vector("list", length = ceiling(n_var / chunk_size))
col_ptr <- 1L

for (start in seq(1L, n_var, by = chunk_size)) {
  end <- min(start + chunk_size - 1L, n_var)
  vidx <- start:end

  seqSetFilter(gds, variant.id = variant_ids[vidx], sample.id = keep_samples, verbose = FALSE)

  gt <- seqGetData(gds, "genotype") # [ploidy x sample x variant]
  if (length(dim(gt)) != 3) stop("Unexpected genotype dimensions.")
  # use first allele for haploid
  gt_mat <- gt[1, , , drop = FALSE]
  gt_mat <- matrix(gt_mat, nrow = length(keep_samples), ncol = length(vidx), byrow = FALSE)

  # mask invalid to NA
  gt_mat[gt_mat < 0] <- NA_integer_

  # local chrom/pos
  chrom_chunk <- chrom[vidx]
  pos_chunk   <- pos[vidx]

  # for each pair, find mismatches where both calls present
  out_list <- vector("list", n_pairs)
  for (p in seq_len(n_pairs)) {
    i <- match(pair_mat[p,1], keep_samples)
    j <- match(pair_mat[p,2], keep_samples)

    gi <- gt_mat[i, ]
    gj <- gt_mat[j, ]
    keep <- !is.na(gi) & !is.na(gj)
    if (!any(keep)) next

    mism <- keep & (gi != gj)
    if (!any(mism)) next

    out_list[[p]] <- data.frame(
      pair = pair_labels[p],
      chromosome = chrom_chunk[mism],
      position   = pos_chunk[mism],
      stringsAsFactors = FALSE
    )
  }

  collect_pairs[[col_ptr]] <- bind_rows(out_list)
  col_ptr <- col_ptr + 1L

  message(sprintf("Chunk %d:%d processed", start, end))
}

seqSetFilter(gds, reset = TRUE)
seqClose(gds)

mismatch_df <- bind_rows(collect_pairs)
if (nrow(mismatch_df) == 0) {
  stop("No between-sample differences found after filtering.")
}

# order factors for clean plotting
mismatch_df$pair <- factor(mismatch_df$pair, levels = sort(unique(mismatch_df$pair)))
mismatch_df$chromosome <- factor(mismatch_df$chromosome, levels = chrom_levels)

# ---- plot: one figure, y = pair, x = position, dot = mismatch ----
# PDF height scales with number of pairs (min 4, max 60 inches)
pdf_height <- max(4, min(60, length(levels(mismatch_df$pair)) * 0.15))

p <- ggplot(mismatch_df, aes(x = position, y = pair)) +
  geom_point(size = 0.15, alpha = 0.35) +
  facet_wrap(~ chromosome, scales = "free_x", ncol = 1) +
  labs(
    title = "Pairwise mitochondrial differences (haploid): dot = sites where samples differ",
    x = "Genomic position",
    y = "Sample pair"
  ) +
  theme_bw(base_size = 9) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey92"),
    axis.text.y = element_text(size = 6)
  )

ggsave(out_pdf, p, width = 12, height = pdf_height, limitsize = FALSE)
message(sprintf("Wrote figure: %s", out_pdf))


























gds_file <- "/scratch/rjp5nc/UK2022_2024/allsites_mito/usdobtusa_mito_allsites.haploid.gds"
out_dir  <- "/scratch/rjp5nc/UK2022_2024/allsites_mito"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- Open GDS ---
g <- seqOpen(gds_file)
samples <- seqGetData(g, "sample.id")
pos     <- seqGetData(g, "position")
n_samp  <- length(samples)
n_var   <- length(pos)
message(sprintf("Loaded %s | %d samples | %d sites", basename(gds_file), n_samp, n_var))

# --- Helper: per-sample depth ---
get_depth_matrix <- function(g, n_var, n_samp) {
  has_dp <- gdsfmt::index.gdsn(g, "annotation/format/DP", silent = TRUE)
  if (!is.null(has_dp)) {
    dp <- seqGetData(g, "annotation/format/DP")
    if (is.list(dp) && !is.null(dp$data)) dp <- dp$data
    dp <- as.matrix(dp)
    if (nrow(dp) == n_var && ncol(dp) == n_samp) dp <- t(dp)
    return(dp)
  }
  has_ad <- gdsfmt::index.gdsn(g, "annotation/format/AD", silent = TRUE)
  if (!is.null(has_ad)) {
    ad <- seqGetData(g, "annotation/format/AD")  # [allele x sample x variant]
    dp <- apply(ad, c(2,3), function(x) sum(x, na.rm = TRUE))
    return(dp)
  }
  stop("No per-sample depth found (neither FORMAT/DP nor FORMAT/AD).")
}

# --- Filter samples by mean DP >= 30 ---
dp <- get_depth_matrix(g, n_var, n_samp)
mean_dp <- rowMeans(dp, na.rm = TRUE)
keep_idx <- which(mean_dp >= 20)
keep_samples <- samples[keep_idx]
k <- length(keep_samples)
if (!k) stop("No samples meet mean depth ≥ 20.")
message(sprintf("Keeping %d/%d samples (mean DP ≥ 20)", k, n_samp))


depth_df <- data.frame(
  sample = samples,
  mean_depth = mean_dp
) %>%
  arrange(desc(mean_depth))

subset(depth_df, sample == "Gilmer5_H5")

# --- Initialize matrices ---
M_mismatch <- matrix(0, k, k)
C_callable <- matrix(0, k, k)

variant_ids <- seq_len(n_var)
chunk_size <- 5000L

# --- Chunked iteration ---
for (start in seq(1L, n_var, by = chunk_size)) {
  end <- min(start + chunk_size - 1L, n_var)
  vids <- variant_ids[start:end]
  seqSetFilter(g, variant.id = vids, sample.id = keep_samples, verbose = FALSE)

  gt <- seqGetData(g, "genotype")  # [ploidy x sample x variant]
  gt <- gt[1, , , drop = FALSE]
  gt <- matrix(gt, nrow = k, ncol = length(vids))
  gt[gt < 0] <- NA_integer_

  G1 <- (gt == 1L) & !is.na(gt)
  G0 <- (gt == 0L) & !is.na(gt)
  U  <- !is.na(gt)
  G1 <- G1 * 1; G0 <- G0 * 1; U <- U * 1

  A <- G1 %*% t(G0)
  M_mismatch <- M_mismatch + A + t(A)
  C_callable <- C_callable + (U %*% t(U))

  message(sprintf("Processed %d–%d of %d variants", start, end, n_var))
}

seqSetFilter(g, reset = TRUE)
seqClose(g)

diag(M_mismatch) <- NA
diag(C_callable) <- NA
prop_diff <- M_mismatch / C_callable


# Symmetrize and handle NA for clustering stability
prop_sym <- (prop_diff + t(prop_diff)) / 2
diag(prop_sym) <- 0
prop_sym[is.na(prop_sym)] <- max(prop_sym, na.rm = TRUE)

# Hierarchical clustering on the symmetrized matrix
hc <- hclust(as.dist(prop_sym), method = "average")
ord <- hc$order
keep_samples_ord <- keep_samples[ord]

# Reorder matrix to clustered order
prop_diff_ord <- prop_diff[ord, ord]

# --- build long df and, crucially, FIX factor ordering for axes ---
heat_df <- as.data.frame(prop_diff_ord)
colnames(heat_df) <- keep_samples_ord
heat_df$sample_i <- keep_samples_ord

library(tidyr)
heat_long <- pivot_longer(
  heat_df, -sample_i, names_to = "sample_j", values_to = "prop_diff"
)

# lock both axes to the clustered order
heat_long$sample_i <- factor(heat_long$sample_i, levels = keep_samples_ord)
heat_long$sample_j <- factor(heat_long$sample_j, levels = keep_samples_ord)

library(ggplot2)
p <- ggplot(heat_long, aes(x = sample_j, y = sample_i, fill = prop_diff)) +
  geom_tile() +
  scale_fill_viridis_c(option = "C", na.value = "white",
                       name = "Proportion\ndifferent",
                       limits = c(0, max(heat_long$prop_diff, na.rm = TRUE))) +
  labs(
    title = "Pairwise mitochondrial differences (haploid)",
    subtitle = "Samples ordered by hierarchical clustering",
    x = NULL, y = NULL
  ) +
  coord_fixed() +
  theme_bw(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
    axis.text.y = element_text(size = 6),
    panel.grid = element_blank(),
    legend.position = "right"
  )

png_file <- file.path(out_dir, "pairwise_mito_diff_heatmap_clustered.png")
pdf_file <- file.path(out_dir, "pairwise_mito_diff_heatmap_clustered.pdf")
ggsave(png_file, p, width = max(8, min(16, length(keep_samples_ord) * 0.18)),
       height = max(8, min(16, length(keep_samples_ord) * 0.18)), dpi = 300)
ggsave(pdf_file, p, width = max(8, min(16, length(keep_samples_ord) * 0.18)),
       height = max(8, min(16, length(keep_samples_ord) * 0.18)))




p <- ggplot(heat_long, aes(x = sample_j, y = sample_i, fill = prop_diff)) +
  geom_tile() +
  scale_fill_gradient(
    low = "white", high = "red",
    na.value = "white",
    name = "Proportion\ndifferent",
    limits = c(0, max(heat_long$prop_diff, na.rm = TRUE))
  ) +
  labs(
    title = "Pairwise mitochondrial differences (haploid)",
    subtitle = "Samples ordered by hierarchical clustering",
    x = NULL, y = NULL
  ) +
  coord_fixed() +
  theme_bw(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
    axis.text.y = element_text(size = 6),
    panel.grid = element_blank(),
    legend.position = "right"
  )

png_file <- file.path(out_dir, "pairwise_mito_diff_heatmap_clustered_20_red.png")
pdf_file <- file.path(out_dir, "pairwise_mito_diff_heatmap_clustered_20_red.pdf")
ggsave(png_file, p, width = max(8, min(16, length(keep_samples_ord) * 0.18)),
       height = max(8, min(16, length(keep_samples_ord) * 0.18)), dpi = 300)
ggsave(pdf_file, p, width = max(8, min(16, length(keep_samples_ord) * 0.18)),
       height = max(8, min(16, length(keep_samples_ord) * 0.18)))