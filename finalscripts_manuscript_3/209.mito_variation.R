#ijob -A berglandlab -c10 -p standard --mem=40G
#module load gcc/11.4.0  openmpi/4.1.4 icu R/4.3.1
#R

#14642bp
library(SeqArray)
library(gdsfmt)
library(dplyr)
library(readr)
library(stringr)
library(purrr)
library(tidyr)
library(ggplot2)

gds_path <- "/scratch/rjp5nc/UK2022_2024/allsites_mito/usdobtusa_mito_allsites.annot2.ALL.gds"
mitotype_csv <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/mito_types.csv"
out_dir <- "/scratch/rjp5nc/UK2022_2024/allsites_mito/unique_snps_by_group"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Mitotype map ----
mitotypes <- read.csv(mitotype_csv, stringsAsFactors = FALSE)
# Normalize columns and drop the first index column if present
stopifnot(all(c("CloneA","Group") %in% names(mitotypes)))
mito <- mitotypes %>%
  transmute(sample = as.character(CloneA),
            group  = as.character(Group)) %>%
  mutate(sample = str_trim(sample), group = str_trim(group)) %>%
  filter(nzchar(sample), nzchar(group))

# ---- Open GDS ----
options(mc.cores = 1)
Sys.setenv(OMP_NUM_THREADS = "1")

g <- seqOpen(gds_path)

## 1) Recompute per-sample mean DP (orientation-safe) and build clean keep_samples
seqResetFilter(g)
DP_all <- seqGetData(g, "annotation/format/DP"); if (is.list(DP_all)) DP_all <- do.call(cbind, DP_all)
g_samples <- seqGetData(g, "sample.id")

if (ncol(DP_all) == length(g_samples)) {
  mean_dp <- colMeans(DP_all, na.rm = TRUE)
} else if (nrow(DP_all) == length(g_samples)) {
  mean_dp <- rowMeans(DP_all, na.rm = TRUE)
} else stop("Unexpected DP matrix shape.")

mask <- (mean_dp > 20) & !is.na(mean_dp)
keep_samples <- g_samples[mask]
keep_samples <- keep_samples[!is.na(keep_samples)]  # drop any NA slots
stopifnot(all(!is.na(keep_samples)))

message(sprintf("Keeping %d/%d (mean DP>20). NA mean-DP samples: %d",
                length(keep_samples), length(g_samples), sum(is.na(mean_dp))))

## 2) Apply filter and lock the *actual* sample order
seqResetFilter(g)
seqSetFilter(g, sample.id = keep_samples, action = "set")
samps <- seqGetData(g, "sample.id")   # definitive order for columns

## 3) Rebuild mito_keep strictly from current filtered order
mitotypes <- read.csv("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/mito_types.csv",
                      stringsAsFactors = FALSE)
new_samples <- c("SRR5012393", "SRR5012394", "SRR5012771",
                 "SRR5012396", "SRR5012770", "SRR5012773")

# Build a data frame to append
new_rows <- data.frame(
  X = seq(max(mitotypes$X, na.rm = TRUE) + 1,
          length.out = length(new_samples)),
  CloneA = new_samples,
  Group = "XX",
  stringsAsFactors = FALSE
)

# Append to mitotypes
mitotypes <- rbind(mitotypes, new_rows)

# Check the new tail
tail(mitotypes)

mitotypes <- subset(mitotypes, Group %in% c("A", "B", "C", "D", "XX"))

mito_keep <- mitotypes %>%
  transmute(sample = trimws(as.character(CloneA)),
            group  = trimws(as.character(Group))) %>%
  semi_join(tibble(sample = samps), by = "sample")

## 4) Pull ALL sites genotypes; ensure GT is [variants x samples]
GT_raw <- seqGetData(g, "genotype")  # haploid: ploidy x variants x samples
if (length(dim(GT_raw)) == 3L) {
  GT <- matrix(GT_raw[1, , ], nrow = dim(GT_raw)[2], ncol = dim(GT_raw)[3])
} else if (is.matrix(GT_raw)) {
  GT <- GT_raw
} else stop("Unexpected GT shape.")
if (ncol(GT) != length(samps) && nrow(GT) == length(samps)) GT <- t(GT)
stopifnot(ncol(GT) == length(samps))

is_missing <- is.na(GT) | (GT >= 3L)

## 5) Choose Group-A reference (highest mean DP among kept A samples)
DP_kept <- seqGetData(g, "annotation/format/DP"); if (is.list(DP_kept)) DP_kept <- do.call(cbind, DP_kept)
if (ncol(DP_kept) == length(samps)) {
  mean_dp_kept <- colMeans(DP_kept, na.rm = TRUE)
} else {
  mean_dp_kept <- rowMeans(DP_kept, na.rm = TRUE)
}
A_samps <- mito_keep %>% filter(group == "A") %>% pull(sample)
stopifnot(length(A_samps) > 0)
ref_sample <- tibble(sample = samps, mean_dp = as.numeric(mean_dp_kept)) %>%
  filter(sample %in% A_samps) %>%
  arrange(desc(mean_dp)) %>% slice(1) %>% pull(sample)
message("Reference sample (Group A): ", ref_sample)

## 6) Count differences per sample vs reference across ALL sites
ref_idx <- match(ref_sample, samps); stopifnot(!is.na(ref_idx))
ref_gt  <- GT[, ref_idx]
ref_mis <- is_missing[, ref_idx]

nsamp <- length(samps)
diff_mat <- matrix(0L, nrow = nsamp, ncol = 2,
                   dimnames = list(samps, c("diffs","compared")))
for (j in seq_len(nsamp)) {
  ok <- !(is_missing[, j] | ref_mis)  # both called
  if (any(ok)) {
    diff_mat[j, "diffs"]    <- sum(GT[ok, j] != ref_gt[ok])
    diff_mat[j, "compared"] <- sum(ok)
  }
}

diff_df <- tibble(
  sample = samps,
  diffs = as.integer(diff_mat[, "diffs"]),
  compared_sites = as.integer(diff_mat[, "compared"]),
  prop_diff = ifelse(compared_sites > 0, diffs / compared_sites, NA_real_)
) %>%
  left_join(mito_keep, by = "sample") %>%
  rename(group = group) %>%
  mutate(is_ref = sample == ref_sample)

## 7) Save + quick plot
out_dir <- "/scratch/rjp5nc/UK2022_2024/allsites_mito/unique_snps_by_group"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_tbl <- file.path(out_dir, paste0("allsite_pairwise_diff_vs_", make.names(ref_sample), ".tsv"))
readr::write_tsv(diff_df, out_tbl)
message("Wrote: ", out_tbl)

p_box <- ggplot(diff_df %>% filter(sample != ref_sample, !is.na(group)),
                aes(x = group, y = prop_diff)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, height = 0, alpha = 0.6) +
  labs(title = paste0("All-site differences vs ", ref_sample, " (Group A ref)"),
       x = "Group", y = "Proportion different (all sites)") +
  theme_bw()
ggsave(file.path(out_dir, paste0("allsite_prop_diff_vs_", make.names(ref_sample), "_by_group.png")),
       p_box, width = 8, height = 5, dpi = 300)




















## ==== 6.5) Prep variant coordinates ====
pos <- as.integer(seqGetData(g, "position"))
chr <- as.character(seqGetData(g, "chromosome"))
stopifnot(length(pos) == nrow(GT), length(chr) == nrow(GT))

# Order chromosomes nicely (optional)
chr_levels <- unique(chr)
chr <- factor(chr, levels = chr_levels)

## ==== 6.6) Windowed differences vs reference ====
window_size <- 200  # 150 ; change as you like
bin <- floor((pos - 1L) / window_size)         # 0-based bin index per site
bin_start <- bin * window_size + 1L
bin_end   <- (bin + 1L) * window_size
bin_mid   <- (bin_start + bin_end) / 2

nsamp <- length(samps)
window_list <- vector("list", nsamp)

for (j in seq_len(nsamp)) {
  ok <- !(is_missing[, j] | ref_mis)              # callable in both
  dif <- as.integer(GT[, j] != ref_gt)            # site-level difference
  dif_ok <- ifelse(ok, dif, 0L)

  df_j <- tibble(
    sample = samps[j],
    chr = chr,
    bin = bin,
    mid = bin_mid,
    ok = as.integer(ok),
    dif_ok = dif_ok
  ) %>%
    group_by(chr, bin, mid, sample) %>%
    summarise(
      compared = sum(ok, na.rm = TRUE),
      diffs    = sum(dif_ok, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(prop_diff = ifelse(compared > 0, diffs / compared, NA_real_))

  window_list[[j]] <- df_j
}

win_df <- bind_rows(window_list) %>%
  left_join(mito_keep, by = c("sample" = "sample")) %>%
  rename(group = group)

# Tidy chromosome for plotting order
win_df <- win_df %>%
  mutate(chr = factor(as.character(chr), levels = as.character(chr_levels)))

## ==== 6.7) Optional: group-weighted mean per window ====
group_mean <- win_df %>%
  group_by(chr, bin, mid, group) %>%
  summarise(
    compared_total = sum(compared, na.rm = TRUE),
    diffs_total    = sum(diffs, na.rm = TRUE),
    prop_mean      = ifelse(compared_total > 0, diffs_total / compared_total, NA_real_),
    .groups = "drop"
  )

## ==== 7) Save tables ====
out_dir <- "/scratch/rjp5nc/UK2022_2024/allsites_mito/unique_snps_by_group"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

readr::write_tsv(win_df, file.path(out_dir, "allsite_windowed_prop_diff_by_sample.tsv"))
readr::write_tsv(group_mean, file.path(out_dir, "allsite_windowed_prop_diff_groupmean.tsv"))

## ==== 8) Plots ====

# (A) Per-sample tracks (thin lines), faceted by chromosome
p_tracks <- ggplot(subset(win_df, !is.na(group)), aes(x = mid, y = prop_diff, group = sample, color = group)) +
  geom_line(alpha = 0.35, linewidth = 0.3, na.rm = TRUE) +
  facet_wrap(~ chr, scales = "free_x", ncol = 1, strip.position = "left") +
  labs(
    title = paste0("Across-genome differences vs ", ref_sample, " (", window_size/1000, " kb windows)"),
    x = "Genomic position (bp)",
    y = "Proportion different",
    color = "Group"
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.spacing.y = unit(0.4, "lines"),
    strip.background = element_rect(fill = "grey90"),
    strip.placement = "outside"
  )

ggsave(file.path(out_dir, paste0("tracks_across_genome_vs_", make.names(ref_sample), "_", window_size/1000, "kb.png")),
       p_tracks, width = 10, height = max(4, length(chr_levels) * 1.2), dpi = 300)

# (B) Cleaner summary: group-weighted mean per window
p_group <- ggplot(subset(group_mean, !is.na(group)), aes(x = mid, y = prop_mean, color = group)) +
  geom_line(linewidth = 0.8, na.rm = TRUE) +
  facet_wrap(~ chr, scales = "free_x", ncol = 1, strip.position = "left") +
  labs(
    title = paste0("Across-genome group mean differences vs ", ref_sample, " (", window_size/1000, " kb windows)"),
    x = "Genomic position (bp)",
    y = "Proportion different (group mean)",
    color = "Group"
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.spacing.y = unit(0.4, "lines"),
    strip.background = element_rect(fill = "grey90"),
    strip.placement = "outside"
  )

ggsave(file.path(out_dir, paste0("groupmean_across_genome_vs_", make.names(ref_sample), "_", window_size/1000, "kb.png")),
       p_group, width = 10, height = max(4, length(chr_levels) * 1.2), dpi = 300)


hist <- ggplot(win_df, aes(x = compared)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "black") +
  labs(
    title = "Distribution of callable sites per window",
    x = "Callable sites per window",
    y = "Count of windows"
  ) +
  theme_bw(base_size = 12)



ggsave(file.path(out_dir, "hist.png"),
       hist, width = 10, height = 4, dpi = 300)












library(SeqArray)
library(gdsfmt)
library(dplyr)
library(tibble)
library(ape)
library(ggtree)
library(ggplot2)

g <- seqOpen(gds_path)
on.exit(seqClose(g), add = TRUE)

## Keep good-coverage samples (reuse your >20 or >30 rule)
seqResetFilter(g)
DP <- seqGetData(g, "annotation/format/DP"); if (is.list(DP)) DP <- do.call(cbind, DP)
samps_all <- seqGetData(g, "sample.id")
mean_dp <- if (ncol(DP) == length(samps_all)) colMeans(DP, na.rm = TRUE) else rowMeans(DP, na.rm = TRUE)
keep <- !is.na(mean_dp) & mean_dp > 20
seqResetFilter(g); seqSetFilter(g, sample.id = samps_all[keep], action = "set")
samps <- seqGetData(g, "sample.id")

## Haploid genotypes: variants x samples
GT_raw <- seqGetData(g, "genotype")
if (length(dim(GT_raw)) == 3L) {
  GT <- matrix(GT_raw[1, , ], nrow = dim(GT_raw)[2], ncol = dim(GT_raw)[3])
} else if (is.matrix(GT_raw)) {
  GT <- GT_raw
} else stop("Unexpected GT shape.")
if (ncol(GT) != length(samps) && nrow(GT) == length(samps)) GT <- t(GT)
stopifnot(ncol(GT) == length(samps))
is_missing <- is.na(GT) | (GT >= 3L)

## ---------- IBS (identity-by-state) ----------
## For haploid 0/1 genotypes: IBS = matches / jointly-callable
n <- ncol(GT)
IBS <- matrix(NA_real_, n, n, dimnames = list(samps, samps))
for (i in seq_len(n)) {
  IBS[i, i] <- 1
  for (j in (i + 1L):n) {
    ok <- !(is_missing[, i] | is_missing[, j])
    denom <- sum(ok)
    if (denom == 0) {
      v <- NA_real_
    } else {
      v <- sum(GT[ok, i] == GT[ok, j]) / denom
    }
    IBS[i, j] <- IBS[j, i] <- v
  }
}
## Convert to a distance for NJ
D <- 1 - IBS

## Replace any NA distances with a large-ish value to avoid nj() failure
if (anyNA(D)) {
  fill <- max(D, na.rm = TRUE)
  D[is.na(D)] <- fill
}

tree <- nj(as.dist(D))

## ---------- Tip metadata: only mito_type (Group), include NAs ----------
mitotypes <- read.csv(mitotype_csv, stringsAsFactors = FALSE)
# (Optional) append SRR... as XX if you want them present with a label
extra <- c("SRR5012393","SRR5012394","SRR5012771","SRR5012396","SRR5012770","SRR5012773")
mitotypes <- bind_rows(
  mitotypes %>% select(CloneA, Group),
  tibble(CloneA = extra, Group = "XX")
)

tip_meta <- tibble(sample = samps) %>%
  left_join(mitotypes %>% transmute(sample = as.character(CloneA),
                                    mito_type = as.character(Group)),
            by = "sample")

tip_meta2 <- tip_meta %>%
  mutate(
    tip = sample,  # rename instead of 'label' to avoid clash
    mito_plot = ifelse(is.na(mito_type), "Unknown", mito_type)
  )

## Verify match with tree tip labels
stopifnot(all(tree$tip.label %in% tip_meta2$tip))

## Palette including 'Unknown'
mt_levels <- sort(unique(tip_meta2$mito_plot))
pal <- setNames(hue_pal()(length(mt_levels)), mt_levels)

## Rectangular tree
p_rect <- ggtree(tree) %<+% tip_meta2 +
  geom_tiplab(aes(label = tip, color = mito_plot), size = 2) +
  scale_color_manual(values = pal, name = "Mitotype") +
  labs(title = "Neighbor-Joining (IBS distance) — mitochondrial all-sites") +
  theme_tree2()

ggsave(file.path(out_dir, "NJ_IBS_rect.png"), p_rect, width = 10, height = 10, dpi = 600)

## Circular tree
p_circ <- ggtree(tree, layout = "circular") %<+% tip_meta2 +
  geom_tiplab(aes(label = tip, color = mito_plot), size = 2) +
  scale_color_manual(values = pal, name = "Mitotype") +
  theme_tree2()

ggsave(file.path(out_dir, "NJ_IBS_circular.png"), p_circ, width = 10, height = 10, dpi = 600)


tree_rooted <- root(tree, outgroup = "SRR5012393", resolve.root = TRUE)

# Optional: drop the outgroup tip from display if you don’t want it shown
# tree_rooted <- drop.tip(tree_rooted, "SRR5012393")

# Then re-run the plotting code on tree_rooted
p_rect <- ggtree(tree_rooted) %<+% tip_meta2 +
  geom_tiplab(aes(label = tip, color = mito_plot), size = 2) +
  scale_color_manual(values = pal, name = "Mitotype") +
  labs(title = "Neighbor-Joining (IBS distance, rooted by SRR5012393)") +
  theme_tree2()

ggsave(file.path(out_dir, "NJ_IBS_rect_rooted.png"), p_rect, width = 10, height = 10, dpi = 600)

p_circ <- ggtree(tree_rooted, layout = "circular") %<+% tip_meta2 +
  geom_tiplab(aes(label = tip, color = mito_plot), size = 2) +
  scale_color_manual(values = pal, name = "Mitotype") +
  labs(title = "Neighbor-Joining (IBS distance, rooted by SRR5012393)") +
  theme_tree2()

ggsave(file.path(out_dir, "NJ_IBS_circular_rooted.png"), p_circ, width = 10, height = 10, dpi = 600)
