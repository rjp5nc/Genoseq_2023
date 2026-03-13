#ijob -A berglandlab -c10 -p standard --mem=40G
#module load gcc/11.4.0  openmpi/4.1.4 icu R/4.3.1; R

#14642bp
# Load required packages
library(SeqArray)
library(SNPRelate)
library(ape)
library(ggtree)
library(ggplot2)
library(gdsfmt)
library(dplyr)
library(tibble)
library(data.table)
library(igraph)
library(SeqVarTools)
library(readr)
library(stringr)
library(purrr)
library(tidyr)
library(ggnewscale)

library(patchwork)
library(foreach)
library(lubridate)




# ---- Step 1: Open GDS file ----



metadata <- read.csv("/project/berglandlab/Robert/UKSequencing2022_2024/old_stuff/2022_2024seqmetadata20250131.csv", header = TRUE)
metadata_with_clone <- read.csv("/project/berglandlab/Robert/UKSequencing2022_2024/old_stuff/metadata_with_clone.csv", header = TRUE)

samples <- read.csv("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/Sample_ID_Species_merged_20251227.csv")

samples <- subset(samples, !grepl("^bdw", Sample_ID_old))

subset(samples, grepl("^bdw", Sample_ID_old))

usobtusasamps <- subset(samples, Species == "Daphnia obtusa")

metadata_with_clone <- read.csv("/project/berglandlab/Robert/UKSequencing2022_2024/old_stuff/metadata_with_clone.csv", header = TRUE)
metadata_with_clone$clone <- trimws(metadata_with_clone$clone)
metadata_with_clone <- subset(metadata_with_clone, clone !="Blank")
metadata_with_clone <- subset(metadata_with_clone, clone !="BLANK")

diffs <- read.csv("/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/mito_types_diff_df.csv")


head(metadata)
metadata$clone <- trimws(metadata$clone)
metadata <- metadata %>% 
  filter(!tolower(clone) %in% c("blank", "blanks", "na", "missing"))

mitotype_csv <- "/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/mito_types.csv"


meta <- fread(mitotype_csv, header = TRUE)
setnames(meta, c("sample.id", "mitotype"))

meta[, sample.id := trimws(as.character(sample.id))]
meta[, mitotype  := trimws(as.character(mitotype))]

meta[tolower(mitotype) %in% c("na","nan","none","unknown",""), mitotype := NA]
meta <- meta[!is.na(mitotype) & sample.id != ""]

# enforce 1:1 mapping
bad <- meta[, .(n = uniqueN(mitotype)), by = sample.id][n > 1]
if (nrow(bad)) stop("ERROR: conflicting mitotypes in CSV")

meta <- meta[, .(mitotype = mitotype[1]), by = sample.id]
meta_samples <- meta$sample.id


seqClose(genofile)
#gds.fn <- "/scratch/rjp5nc/UK2022_2024/mito_vcf/usdobtusa_mito_genotyped.gds"             # output GDS file
gds.fn <- "/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/usdobtusa_mito_allsites_all.haploid.annotated2.gds"
genofile <- seqOpen(gds.fn)

seqResetFilter(genofile)



dp_filt <- seqGetData(genofile, "annotation/format/DP")  # matrix: variants x samples
sids_filt <- seqGetData(genofile, "sample.id")

# Mean depth per sample across filtered variants
mean_dp_per_sample <- rowMeans(dp_filt, na.rm = TRUE)

# (Optional) also report how many variants contributed (non-missing DP) per sample
n_dp_nonmissing <- rowSums(!is.na(dp_filt))

# Build output table
out_df <- data.frame(
  sample = sids_filt,
  mean_depth = mean_dp_per_sample,
  n_sites_nonmissing = n_dp_nonmissing,
  stringsAsFactors = FALSE
)

# Write file
out_path <- "/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/sample_mean_depth.csv"
write.csv(out_df, out_path, row.names = FALSE)

out_path

library(dplyr)
library(readr)

# out_df already exists with columns: sample, mean_depth, n_sites_nonmissing

groups <- read_csv("/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/mito_types.csv") %>%
  rename(sample = sampleA)

top_per_group <- out_df %>%
  inner_join(groups, by = "sample") %>%
  group_by(Group) %>%
  slice_max(order_by = mean_depth, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(Group)

top_per_group
write_csv(top_per_group,
          "/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/top_sample_per_group_by_DP.csv")

# also write just the samples for bash
writeLines(top_per_group$sample,
           "/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/top_samples.txt")





samples_to_keep <- usobtusasamps %>% filter(Species == "Daphnia obtusa") %>% pull(Sample_ID)
#samples_to_keep <- metadata_with_clone %>% filter(accuratelocation == "P63") %>% pull(Well)
#samples_to_keep <- metadata_with_clone %>% filter(location == "UK" | accuratelocation == "P759") %>% pull(Well)

#samples_to_keep <- metadata_with_clone$Well

unique(metadata_with_clone$accuratelocation)

seqSetFilter(genofile, sample.id = samples_to_keep)

# ---- Step 2: Filter variants with missing rate < 0.05 ----
miss_rate_per_sample <- seqMissing(genofile, per.variant = FALSE)
miss_rate_per_variant <- seqMissing(genofile, per.variant = TRUE)
sample_ids <- seqGetData(genofile, "sample.id")
valid_samples <- sample_ids[miss_rate_per_sample < 0.10]
miss_rate_per_variant <- seqMissing(genofile, per.variant=TRUE)
valid_variants <- seqGetData(genofile, "variant.id")[miss_rate_per_variant < 0.10]


final_valid_samples <- intersect(valid_samples, samples_to_keep)

# seqSetFilter(genofile, sample.id = final_valid_samples)

miss_rate <- seqMissing(genofile, per.variant = TRUE)
dp <- seqGetData(genofile, "annotation/format/DP")
mean_depth <- rowMeans(dp, na.rm = TRUE)
keep <- which(miss_rate < 0.15)
# keep <- which(miss_rate < 0.10)


#ADD THIS FOR COI ANALYSIS - sites from Penton 2004, COI primers aligned to Ref

# keep <- keep[
#   keep >= 1421 & keep <= 2129
# ]


g_samples <- seqGetData(genofile, "sample.id")

mask <- (mean_depth > 30) & !is.na(mean_depth) 
keep_samples <- g_samples[mask]
keep_samples <- keep_samples[!is.na(keep_samples)]  # drop any NA slots
stopifnot(all(!is.na(keep_samples)))
final_valid_samples <- intersect(keep_samples, final_valid_samples)

samps300 <- subset(diffs, diffs < 500)
samps300samps <- samps300$sample

final_valid_samples <- intersect(samps300samps, final_valid_samples)

seqResetFilter(genofile)

keep_samples <- unique(c(keep_samples, "RobertUK_B10"))

seqSetFilter(genofile, sample.id = keep_samples, variant.id = keep)


# ------------------------------------------------------------
# 1) Samples
# ------------------------------------------------------------
samps <- seqGetData(genofile, "sample.id")
nsamp <- length(samps)

# ------------------------------------------------------------
# 2) Genotypes (variants x samples)
# ------------------------------------------------------------
GT_raw <- seqGetData(genofile, "genotype")  # haploid: ploidy x sample x variant  (usually)

stopifnot(length(dim(GT_raw)) == 3L)

d <- dim(GT_raw)
# d[1]=ploidy; d[2] and d[3] are sample/variant but can be swapped depending on file

# We know how many samples we asked for:
nsamp <- length(seqGetData(genofile, "sample.id"))

# If second dim matches nsamp, it's ploidy x samples x variants
if (d[2] == nsamp) {
  GT <- t(GT_raw[1, , ])        # (variants x samples)
# If third dim matches nsamp, it's ploidy x variants x samples
} else if (d[3] == nsamp) {
  GT <- GT_raw[1, , ]           # already (variants x samples)
} else {
  stop("Cannot infer genotype dimension order: dim(genotype) = ",
       paste(d, collapse=" x "), " ; nsamp = ", nsamp)
}

GT <- as.matrix(GT)
colnames(GT) <- seqGetData(genofile, "sample.id")
nvar <- nrow(GT)

# Missing genotypes
is_missing <- is.na(GT) | (GT >= 3L)

message("GT dim = ", paste(dim(GT), collapse=" x "), " (variants x samples)")












## 5) Choose Group-A reference (highest mean DP among kept A samples)
DP_kept <- seqGetData(genofile, "annotation/format/DP"); if (is.list(DP_kept)) DP_kept <- do.call(cbind, DP_kept)
if (ncol(DP_kept) == length(samps)) {
  mean_dp_kept <- colMeans(DP_kept, na.rm = TRUE)
} else {
  mean_dp_kept <- rowMeans(DP_kept, na.rm = TRUE)
}



ref_sample <- "Gilmer5_B4"


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

diff_mat2 <- as.data.frame(diff_mat)

diff_df <- tibble(
  sample = samps,
  diffs = as.integer(diff_mat[, "diffs"]),
  compared_sites = as.integer(diff_mat[, "compared"]),
  prop_diff = ifelse(compared_sites > 0, diffs / compared_sites, NA_real_)
) %>%
  left_join(meta, by = c("sample" = "sample.id")) %>%   # <-- fix join
  rename(group = mitotype) %>%                          # <-- mitotype -> group
  mutate(is_ref = sample == ref_sample)

## 7) Save + quick plot
out_dir <- "/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/unique_snps_by_group"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_tbl <- file.path(out_dir, paste0("allsite_pairwise_diff_vs_", make.names(ref_sample), ".tsv"))
readr::write_tsv(diff_df, out_tbl)
message("Wrote: ", out_tbl)



plot_df <- diff_df %>%
  mutate(group = str_trim(as.character(group))) %>%
  filter(
    !is.na(group),
    group != "",
    !tolower(group) %in% c("na", "n/a"),
    sample != ref_sample   # drop the ref sample no matter what
  )

p_box <- ggplot(plot_df, aes(x = group, y = prop_diff)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, height = 0, alpha = 0.6) +
  labs(
    title = paste0("All-site differences vs ", ref_sample, " (Group A ref)"),
    x = "Group", y = "Proportion different (all sites)"
  ) +
  theme_bw()

ggsave(
  file.path(out_dir, paste0("allsite_prop_diff_vs_", make.names(ref_sample), "_by_group.pdf")),
  p_box, width = 6, height = 5, dpi = 600
)




p_box <- ggplot(diff_df %>% filter(sample != ref_sample),
                aes(x = group, y = prop_diff)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, height = 0, alpha = 0.6) +
  labs(title = paste0("All-site differences vs ", ref_sample, " (Group A ref)"),
       x = "Group", y = "Proportion different (all sites)") +
  theme_bw()
ggsave(file.path(out_dir, paste0("allsite_prop_diff_vs_", make.names(ref_sample), "_by_group_WithNAsover20dp.png")),
       p_box, width = 4, height = 5, dpi = 300)

write.csv(diff_df,"/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/mito_types_diff_df.csv")




















# ref_sample <- "SRR14204327"


# ## 6) Count differences per sample vs reference across ALL sites
# ref_idx <- match(ref_sample, samps); stopifnot(!is.na(ref_idx))
# ref_gt  <- GT[, ref_idx]
# ref_mis <- is_missing[, ref_idx]

# nsamp <- length(samps)
# diff_mat <- matrix(0L, nrow = nsamp, ncol = 2,
#                    dimnames = list(samps, c("diffs","compared")))
# for (j in seq_len(nsamp)) {
#   ok <- !(is_missing[, j] | ref_mis)  # both called
#   if (any(ok)) {
#     diff_mat[j, "diffs"]    <- sum(GT[ok, j] != ref_gt[ok])
#     diff_mat[j, "compared"] <- sum(ok)
#   }
# }

# diff_mat2 <- as.data.frame(diff_mat)

# diff_df <- tibble(
#   sample = samps,
#   diffs = as.integer(diff_mat[, "diffs"]),
#   compared_sites = as.integer(diff_mat[, "compared"]),
#   prop_diff = ifelse(compared_sites > 0, diffs / compared_sites, NA_real_)
# ) %>%
#   left_join(meta, by = c("sample" = "sample.id")) %>%   # <-- fix join
#   rename(group = mitotype) %>%                          # <-- mitotype -> group
#   mutate(is_ref = sample == ref_sample)

# ## 7) Save + quick plot
# out_dir <- "/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/unique_snps_by_group"
# dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
# out_tbl <- file.path(out_dir, paste0("allsite_pairwise_diff_vs_", make.names(ref_sample), ".tsv"))
# readr::write_tsv(diff_df, out_tbl)
# message("Wrote: ", out_tbl)

# p_box <- ggplot(diff_df %>% filter(sample != ref_sample, !is.na(group)),
#                 aes(x = group, y = prop_diff)) +
#   geom_boxplot(outlier.shape = NA) +
#   geom_jitter(width = 0.15, height = 0, alpha = 0.6) +
#   labs(title = paste0("All-site differences vs ", ref_sample, " (Group A ref)"),
#        x = "Group", y = "Proportion different (all sites)") +
#   theme_bw()
# ggsave(file.path(out_dir, paste0("allsite_prop_diff_vs_", make.names(ref_sample), "_by_group.png")),
#        p_box, width = 6, height = 5, dpi = 600)

# p_box <- ggplot(diff_df %>% filter(sample != ref_sample),
#                 aes(x = group, y = prop_diff)) +
#   geom_boxplot(outlier.shape = NA) +
#   geom_jitter(width = 0.15, height = 0, alpha = 0.6) +
#   labs(title = paste0("All-site differences vs ", ref_sample, " (Group A ref)"),
#        x = "Group", y = "Proportion different (all sites)") +
#   theme_bw()
# ggsave(file.path(out_dir, paste0("allsite_prop_diff_vs_", make.names(ref_sample), "_by_group_WithNAsover20dp.png")),
#        p_box, width = 4, height = 5, dpi = 300)

# write.csv(diff_df,"/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/mito_types_diff_df_SRR14204327.csv")














# ------------------------------------------------------------
# Pick genetically-unique samples: collapse identical genotypes,
# keep the one with highest mean depth (within current filters)
# ------------------------------------------------------------

# 0) Load digest for fast hashing (recommended)
# install.packages("digest")  # if needed
library(digest)
library(data.table)

# 1) Mean depth per sample on the CURRENT filtered data
DP_raw <- seqGetData(genofile, "annotation/format/DP")
if (is.list(DP_raw)) DP_raw <- do.call(cbind, DP_raw)  # usually variants x samples
DP_mat <- as.matrix(DP_raw)

# Align DP_mat to (variants x samples) like GT
# GT is (variants x samples)
if (nrow(DP_mat) == ncol(GT) && ncol(DP_mat) == nrow(GT)) {
  # DP is (samples x variants) -> transpose
  DP_mat <- t(DP_mat)
} else if (!(nrow(DP_mat) == nrow(GT) && ncol(DP_mat) == ncol(GT))) {
  stop("DP matrix dims do not match GT dims. DP: ",
       paste(dim(DP_mat), collapse="x"),
       " ; GT: ", paste(dim(GT), collapse="x"))
}

mean_dp_by_sample <- colMeans(DP_mat, na.rm = TRUE)
stopifnot(length(mean_dp_by_sample) == ncol(GT))
names(mean_dp_by_sample) <- colnames(GT)

# 2) Build a genotype "signature" per sample
# Replace missing with a sentinel so missing patterns are part of identity
GT_filled <- GT
GT_filled[is_missing] <- -9L

# Hash each sample's genotype vector (fast + avoids huge strings)
geno_hash <- vapply(
  seq_len(ncol(GT_filled)),
  function(j) digest(GT_filled[, j], algo = "xxhash64"),
  character(1)
)
names(geno_hash) <- colnames(GT_filled)

# 3) Group by identical genotype hash
dt <- data.table(
  sample = names(geno_hash),
  hash   = unname(geno_hash),
  meanDP = as.numeric(mean_dp_by_sample[names(geno_hash)])
)

# 4) Choose representative per hash = highest meanDP
# (tie-breaker: keep first alphabetically for determinism)
setorder(dt, hash, -meanDP, sample)

rep_dt <- dt[, .(
  representative = sample[1],
  rep_meanDP     = meanDP[1],
  n_in_group     = .N,
  members        = paste(sample, collapse = ",")
), by = hash]

unique_samples2 <- rep_dt$representative

# 5) Save outputs
out_dir <- "/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/unique_snps_by_group"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fwrite(rep_dt, file.path(out_dir, "genotype_duplicate_groups.csv"), sep = "\t")

fwrite(
  data.table(sample.id = unique_samples2),
  file.path(out_dir, "genetically_unique_samples_highDP.txt"),
  col.names = FALSE
)

message("Unique genetically (kept): ", length(unique_samples2),
        " out of ", ncol(GT), " samples")









## ==== 6.5) Prep variant coordinates ====
pos <- as.integer(seqGetData(genofile, "position"))
chr <- as.character(seqGetData(genofile, "chromosome"))
stopifnot(length(pos) == nrow(GT), length(chr) == nrow(GT))

# Order chromosomes nicely (optional)
chr_levels <- unique(chr)
chr <- factor(chr, levels = chr_levels)

## ==== 6.6) Windowed differences vs reference ====
window_size <- 50  # 150 ; change as you like
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
  left_join(meta, by = c("sample" = "sample.id")) %>%
  rename(group = mitotype)

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













# ------------------------------------------------------------
# 3) Pairwise differences (ALL vs ALL)
# ------------------------------------------------------------
diff_counts <- matrix(0L, nsamp, nsamp, dimnames = list(samps, samps))
compared    <- matrix(0L, nsamp, nsamp, dimnames = list(samps, samps))

for (i in seq_len(nsamp)) {
  for (j in i:nsamp) {
    ok <- !(is_missing[, i] | is_missing[, j])
    n_ok <- sum(ok)

    compared[i, j] <- compared[j, i] <- n_ok

    if (n_ok > 0) {
      d <- sum(GT[ok, i] != GT[ok, j])
      diff_counts[i, j] <- diff_counts[j, i] <- d
    }
  }
}

diag(diff_counts) <- 0L
diag(compared) <- 0L

# ------------------------------------------------------------
# 4) Proportion differences & similarity
# ------------------------------------------------------------
prop_diff <- diff_counts / compared
Similarity <- 1 - prop_diff

# ------------------------------------------------------------
# 5) Long table (SAFE, aligned)
# ------------------------------------------------------------
DT_prop <- as.data.table(as.table(prop_diff))
setnames(DT_prop, c("sampleA","sampleB","prop_diff"))

DT_diff <- as.data.table(as.table(diff_counts))
setnames(DT_diff, c("sampleA","sampleB","diffs"))

DT_comp <- as.data.table(as.table(compared))
setnames(DT_comp, c("sampleA","sampleB","compared_sites"))

pairwise_long <- merge(DT_prop, DT_diff, by = c("sampleA","sampleB"))
pairwise_long <- merge(pairwise_long, DT_comp, by = c("sampleA","sampleB"))

pairwise_long[, Similarity := 1 - prop_diff]
pairwise_long <- pairwise_long[sampleA != sampleB]

pairwise_long[, Similarity := 1 - prop_diff]
pairwise_long <- pairwise_long[sampleA != sampleB]

# ------------------------------------------------------------
# 6) Save pairwise outputs
# ------------------------------------------------------------
out_dir <- "/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/unique_snps_by_group"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fwrite(pairwise_long,
       file.path(out_dir, "allsite_pairwise_long.tsv"))

saveRDS(
  list(diff_counts = diff_counts,
       compared = compared,
       prop_diff = prop_diff),
  file = file.path(out_dir, "allsite_pairwise_matrices.rds")
)





pairwise_long <- as.data.table(pairwise_long)
meta <- as.data.table(meta)

# -----------------------------
# 1) Deduplicate by sample IDs
# -----------------------------
pairwise_uniq <- pairwise_long[sampleA != sampleB]

pairwise_uniq[, sample_min := pmin(sampleA, sampleB)]
pairwise_uniq[, sample_max := pmax(sampleA, sampleB)]

# keep one row per unordered pair
pairwise_uniq <- unique(pairwise_uniq, by = c("sample_min", "sample_max"))

# --------------------------------
# 2) Join mitotypes onto each side
# --------------------------------
pairwise_anno <- merge(
  pairwise_uniq,
  meta[, .(sampleA = sample.id, mitoA = mitotype)],
  by = "sampleA",
  all.x = TRUE
)

pairwise_anno <- merge(
  pairwise_anno,
  meta[, .(sampleB = sample.id, mitoB = mitotype)],
  by = "sampleB",
  all.x = TRUE
)

# ----------------------------------------
# 3) Labels (directed + undirected version)
# ----------------------------------------
pairwise_anno[, mito_comp := paste0(mitoA, "_", mitoB)]
pairwise_anno[, mito_comp_sorted := paste0(pmin(mitoA, mitoB), "_", pmax(mitoA, mitoB))]

# ----------------------------------------
# (Optional) drop NA / "NA" mitotypes
# ----------------------------------------
pairwise_anno <- pairwise_anno[
  !is.na(mitoA) & !is.na(mitoB) & mitoA != "NA" & mitoB != "NA"
]

# sanity check: no duplicate unordered pairs remain
stopifnot(pairwise_anno[, .N, by = .(sample_min, sample_max)][N > 1, .N] == 0)

pairwise_anno[, mito_comp := paste0(pmin(mitoA, mitoB), "_", pmax(mitoA, mitoB))]

table(pairwise_anno$mito_comp)

head(pairwise_anno)






p_box2 <- ggplot(pairwise_anno, aes(x = mito_comp, y = prop_diff)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, height = 0, alpha = 0.6) +
  labs(
    title = paste0("Mitotype comparisons"),
    x = "Group", y = "Proportion different (all sites)"
  ) + xlab("Mitotype Comparison")+
  theme_bw()

ggsave(
  file.path(out_dir, paste0("allsite_prop_diff_all.pdf")),
  p_box2, width = 15, height = 5, dpi = 600
)



p_box3 <- ggplot(subset(pairwise_anno, mito_comp == "A_A" | 
        mito_comp == "A_B" |
         mito_comp == "A_C" |
          mito_comp == "B_B" |          
          mito_comp == "B_C" |
                    mito_comp == "C_C" ), aes(x = mito_comp, y = prop_diff)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, height = 0, alpha = 0.6) +
  labs(
    title = paste0("NA2 to NA2"),
    x = "Group", y = "Proportion different (all sites)"
  ) + ylim(0,0.05)+
  theme_bw()

ggsave(
  file.path(out_dir, paste0("allsite_prop_diff_ABC.pdf")),
  p_box3, width = 6, height = 5, dpi = 600
)




p_box4 <- ggplot(subset(pairwise_anno, mito_comp == "D_D" | 
        mito_comp == "D_E" |
         mito_comp == "D_F" |
          mito_comp == "E_E" |          
          mito_comp == "E_F" |
                    mito_comp == "F_F" ), aes(x = mito_comp, y = prop_diff)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, height = 0, alpha = 0.6) +
  labs(
    title = paste0("NA1 to NA1"),
    x = "Group", y = "Proportion different (all sites)"
  ) + ylim(0,0.05)+
  theme_bw()

ggsave(
  file.path(out_dir, paste0("allsite_prop_diff_DEF.pdf")),
  p_box4, width = 6, height = 5, dpi = 600
)


p_box5 <- p_box4 + p_box3 +
  plot_annotation(tag_levels = "A")


ggsave(
  file.path(out_dir, paste0("allsite_prop_diff_NA1vNA2.pdf")),
  p_box5, width = 10, height = 5, dpi = 600
)

# ------------------------------------------------------------
# 7) Distance matrix for NJ tree
# ------------------------------------------------------------
dist_matrix <- as.dist(prop_diff)

# ------------------------------------------------------------
# 8) Build & root tree
# ------------------------------------------------------------
tree <- nj(dist_matrix)
tree_rooted <- root(tree,
                    outgroup = "RobertUK_B10",
                    resolve.root = TRUE)




pairwise_long_diff_plots <- ggplot(pairwise_long, aes(x = diffs)) +
  geom_histogram(bins = 100) +
  theme_bw() + xlim(0,100)+
  xlab("diffs") +
  ylab("Frequency")



ggsave(file.path(out_dir, "diff_rates.png"),
       pairwise_long_diff_plots, width = 7, height = 6, dpi = 300)




pairwise_long_diff_plots_all <- ggplot(pairwise_long, aes(x = diffs)) +
  geom_histogram(bins = 100) +
  theme_bw() + 
  xlab("diffs") +
  ylab("Frequency")



ggsave(file.path(out_dir, "diff_rates_all.png"),
       pairwise_long_diff_plots_all, width = 7, height = 6, dpi = 300)





pairwise_long_diff_plots_300 <- ggplot(pairwise_long, aes(x = diffs)) +
  geom_histogram(bins = 300) +
  theme_bw() + xlim(0,300)+
  xlab("diffs") +
  ylab("Frequency")



ggsave(file.path(out_dir, "diff_rates_300.png"),
       pairwise_long_diff_plots_300, width = 7, height = 6, dpi = 300)



pairwise_long_diff_plots_300_perc <- ggplot(pairwise_long, aes(x = Similarity)) +
  geom_histogram(bins = 100) +
    geom_vline(xintercept = 0.995, linetype = "dashed") +
  theme_bw() + xlim(0.85,1)+
  xlab("diffs") +
  ylab("Frequency")



ggsave(file.path(out_dir, "diff_rates_300_perc.png"),
       pairwise_long_diff_plots_300_perc, width = 7, height = 6, dpi = 300)


ggsave(file.path(out_dir, "diff_rates_300_perc.pdf"),
       pairwise_long_diff_plots_300_perc, width = 7, height = 6, dpi = 300)




diffs15 <- subset(pairwise_long, Similarity >= 0.995)


diffs15 <- setDT(diffs15)

# 3. Build graph
g <- graph_from_data_frame(diffs15[, c("sampleA", "sampleB")], directed = FALSE)

png("/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/grouping_graph_mito.png",
    width = 2000, height = 2000, res = 300)
plot(g, vertex.size = 5, vertex.label = NA)  # basic igraph plot
dev.off()


# 4. Find connected components
comp <- components(g)

# 5. Assign group letters
n_groups <- max(comp$membership)
group_levels <- sprintf("%s", sapply(1:n_groups, function(i) {
  # Excel-style base-26
  lab <- ""
  x <- i
  while (x > 0) {
    x <- x - 1
    lab <- paste0(LETTERS[(x %% 26) + 1], lab)
    x <- x %/% 26
  }
  lab
}))

group_letters <- setNames(group_levels[comp$membership], names(comp$membership))

# 6. Add group column (both clones exist in the graph)
diffs15_samp <- diffs15 %>%
  mutate(Group = group_letters[sampleA])

# ------------------------------------------------------------
# 9) Group metadata
# ------------------------------------------------------------
unique_clones <- diffs15_samp %>%
  distinct(sampleA, Group)


# ------------------------------------------------------------
# 10) Tip colors (CORRECT alignment)
# ------------------------------------------------------------
groups <- sort(unique(unique_clones$Group))

group_colors <- c(
  A       = "#F8766D",
  B       = "#C49A00",
  C       = "#53B400",
  D       = "#00C094",
  E       = "#00B6EB",
  F       = "#A58AFF",
  Unknown = "#FB61D7"
)


tip_groups <- unique_clones$Group[
  match(tree_rooted$tip.label, unique_clones$sample)
]

tip_colors <- group_colors[tip_groups]
tip_colors[is.na(tip_colors)] <- "grey50"


samples$preferred_label <- ifelse(
  !is.na(samples$Sample_ID_old) & samples$Sample_ID_old != "",
  samples$Sample_ID_old,
  samples$Sample_ID
)

label_map <- setNames(samples$preferred_label, samples$Sample_ID)

tree_rooted$tip.label <- ifelse(
  tree_rooted$tip.label %in% names(label_map),
  label_map[tree_rooted$tip.label],
  tree_rooted$tip.label
)

## ---- Pass 2: if in metadata_with_clone and accuratelocation is present, add it ----
meta_loc <- metadata_with_clone

meta_loc$accuratelocation <- trimws(meta_loc$accuratelocation)
meta_loc$has_loc <- !is.na(meta_loc$accuratelocation) & meta_loc$accuratelocation != ""

loc_map <- setNames(meta_loc$accuratelocation[meta_loc$has_loc],
                    meta_loc$Well[meta_loc$has_loc])

tree_rooted$tip.label <- ifelse(
  tree_rooted$tip.label %in% names(loc_map),
  paste0(loc_map[tree_rooted$tip.label]),
  tree_rooted$tip.label
)

tree_rooted2 <- drop.tip(tree_rooted, "bdw4")
tree_rooted2$edge.length[tree_rooted2$edge.length < 0] <- 0

# ------------------------------------------------------------
# 11) Plot tree
# ------------------------------------------------------------
png("/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/usobtusa_mito_groups.png",
    res = 300, width = 4000, height = 9000)

plot.phylo(tree_rooted2,
           type = "phylogram",
           cex = 0.8,
           label.offset = 0.01,
           no.margin = TRUE,
           tip.color = tip_colors,
           main = "Rooted NJ Tree (mitochondrial, all sites)")

legend("topleft",
       legend = groups,
       col = group_colors[groups],
       pch = 19,
       pt.cex = 1.5,
       cex = 1,
       bty = "n",
       title = "Group")

dev.off()



xr <- range(ape::node.depth.edgelength(tree_rooted2))
xmax <- quantile(ape::node.depth.edgelength(tree_rooted2), 0.98)  # adjust 0.95–0.995





png("/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/usobtusa_mito_groups_circ.png",
    res = 300, width = 15000, height = 9000)

plot.phylo(tree_rooted2,
           type = "fan",
           cex = 0.8,
           no.margin = TRUE,
           tip.color = tip_colors,
           x.lim = c(0, xmax),
           main = "Rooted NJ Tree (mitochondrial, all sites)")

legend("topleft",
       legend = groups,
       col = group_colors[groups],
       pch = 19,
       pt.cex = 1.5,
       cex = 1,
       bty = "n",
       title = "Group")

dev.off()


pdf(
  "/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/usobtusa_mito_groups_circ.pdf",
  width = 20,
  height = 15)

plot.phylo(
  tree_rooted2,
  type = "fan",
  cex = 0.8,
  no.margin = TRUE,
  tip.color = tip_colors,
  x.lim = c(0, xmax),
  main = "Rooted NJ Tree (mitochondrial, all sites)"
)

legend(
  "topleft",
  legend = groups,
  col = group_colors[groups],
  pch = 19,
  pt.cex = 1.5,
  cex = 1,
  bty = "n",
  title = "Group"
)

dev.off()



tree_rooted <- root(tree, outgroup = "RobertUK_B10", resolve.root = TRUE)
tree_id <- tree_rooted
tree_id$edge.length[tree_id$edge.length < 0] <- 0

# Drop bdw4 from the TREE first
tree_id <- drop.tip(tree_id, "SRR5012393")

# Build display labels (now aligned with tree tips)
display_label <- tree_id$tip.label
display_label <- ifelse(display_label %in% names(label_map),
                        label_map[display_label],
                        display_label)
display_label <- ifelse(display_label %in% names(loc_map),
                        loc_map[display_label],
                        display_label)

# Build tip_df (keyed by original sample IDs in tree_id)
prefix_pat <- "^(RAP|bdw|PYR|FS|EBG)"
tip_df <- data.frame(
  label = tree_id$tip.label,
  display_label = display_label,
  mitotype = diffs$group[match(tree_id$tip.label, diffs$sample)],
  tip_class = ifelse(grepl(prefix_pat, display_label, ignore.case = TRUE),
                     "RAP/bdw/PYR/FS/EBG", "Other"),
  stringsAsFactors = FALSE
)

# Optional sanity checks
stopifnot(!("bdw4" %in% tree_id$tip.label))
stopifnot(!("bdw4" %in% tip_df$label))

# Colors for mitotypes
mitos <- sort(unique(na.omit(tip_df$mitotype)))
mito_cols <- setNames(rainbow(length(mitos)), mitos)

# Cap edges for plotting only
xmax <- as.numeric(quantile(ape::node.depth.edgelength(tree_id), 0.98))
tree_plot <- tree_id
tree_plot$edge.length <- pmin(tree_plot$edge.length, xmax / 5)

# Plot
p_circ <- ggtree(tree_plot, layout = "fan") %<+% tip_df +
  geom_tiplab(
    aes(label = display_label, color = mitotype),
    size = 2,
    offset = 0.0015   # <- adjust this value if needed
  ) +
  scale_color_manual(values = group_colors, name = "Mitotype", na.value = "grey70") +
  ggnewscale::new_scale_color() +
  geom_tippoint(aes(color = tip_class), size = 2.5) +
  scale_color_manual(
    values = c("RAP/bdw/PYR/FS/EBG" = "black", "Other" = "grey60"),
    name = "Prefix group"
  ) +
  theme_void()

ggsave(
  "/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/usobtusa_mito_groups_dots_circ.png",
  p_circ, width = 12, height = 12, dpi = 300
)

ggsave(
  "/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/usobtusa_mito_groups_dots_circ.pdf",
  p_circ, width = 12, height = 12, dpi = 300
)

p_rect <- ggtree(tree_rooted2, layout = "rectangular") %<+% tip_df +
  geom_tippoint(aes(color = tip_class), size = 2.5) +
  scale_color_manual(values = c(
    "RAP/bdw/PYR/FS/EBG" = "black",
    "Other" = "grey60"
  )) +
  coord_cartesian(xlim = c(0, xmax), clip = "on") +
  theme_bw() +
  labs(color = "Tip group")


ggsave(
  "/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/usobtusa_mito_groups_dots_rect.png",
  p_rect,
  width = 12,
  height = 12,
  dpi = 300
)



# ------------------------------------------------------------
# 12) Save tree + metadata
# ------------------------------------------------------------
write.tree(tree_rooted2,
           file = "/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/usdobtusa_tree_nwk.nwk")

fwrite(unique_clones,
       "/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/mito_types.csv")

# ------------------------------------------------------------
# 13) Sanity checks
# ------------------------------------------------------------
message("Samples: ", nsamp)
message("Variants: ", nvar)
summary(pairwise_long$compared_sites)

summary(subset(pairwise_long, compared_sites < 10000)$compared_sites)
table(subset(pairwise_long, compared_sites < 10000)$sampleB)

SRR14204909
RobertUK_B10

seqClose(genofile)










suppressPackageStartupMessages({
  library(SeqArray)
  library(data.table)
  library(ggplot2)
})

# ----------------------------
# Paths
# ----------------------------
gds_path     <- "/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/usdobtusa_mito_allsites_all.haploid.annotated2.gds"
mitotype_csv <- "/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/mito_types.csv"
out_dir      <- "/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/unique_snps_by_group"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# Parameters
# ----------------------------
dp_cutoff <- 30
miss_cutoff <- 0.05
maf_cutoff <- 0.00
impute_missing <- TRUE

# ----------------------------
# Read + clean metadata
# KEEP A, DROP ALL UNASSIGNED
# ----------------------------
meta <- fread(mitotype_csv, header = TRUE)
setnames(meta, c("sample.id", "mitotype"))

meta[, sample.id := trimws(as.character(sample.id))]
meta[, mitotype  := trimws(as.character(mitotype))]

meta[tolower(mitotype) %in% c("na","nan","none","unknown",""), mitotype := NA]
meta <- meta[!is.na(mitotype) & sample.id != ""]

# enforce 1:1 mapping
bad <- meta[, .(n = uniqueN(mitotype)), by = sample.id][n > 1]
if (nrow(bad)) stop("ERROR: conflicting mitotypes in CSV")

meta <- meta[, .(mitotype = mitotype[1]), by = sample.id]
meta_samples <- meta$sample.id

# ----------------------------
# Open GDS
# ----------------------------
gds <- seqOpen(gds_path)

# ----------------------------
# Sample filter 1: must have mitotype
# ----------------------------
gds_samples <- seqGetData(gds, "sample.id")
seqSetFilter(gds, sample.id = gds_samples, verbose = TRUE)

# ----------------------------
# Sample filter 2: mean DP > cutoff
# ----------------------------
DP <- seqGetData(gds, "annotation/format/DP")
DP_mat <- if (is.list(DP)) do.call(rbind, DP) else DP

g_samples <- seqGetData(gds, "sample.id")

mean_dp <- if (ncol(DP_mat) == length(g_samples)) {
  colMeans(DP_mat, na.rm = TRUE)
} else {
  rowMeans(DP_mat, na.rm = TRUE)
}

keep_samples <- g_samples[!is.na(mean_dp) & mean_dp > dp_cutoff]







seqSetFilter(gds, sample.id = final_valid_samples, verbose = TRUE)

# LOCK sample IDs (FILTER-AWARE)
sample_ids <- seqGetData(gds, "sample.id")
message("Final samples kept: ", length(sample_ids))

# HARD CHECK
if (anyNA(match(sample_ids, meta$sample.id))) {
  stop("ERROR: sample without mitotype survived filtering")
}

# ----------------------------
# Read genotype (KNOWN LAYOUT: 1 x sample x variant)
# ----------------------------
geno <- seqGetData(gds, "genotype")
stopifnot(dim(geno)[1] == 1, dim(geno)[2] == length(sample_ids))

X <- geno[1, , , drop = TRUE]   # sample x variant
rownames(X) <- sample_ids
X[!is.na(X) & X > 1] <- 1

message(sprintf("Loaded X: %d samples x %d variants", nrow(X), ncol(X)))

# ----------------------------
# Variant filters (in memory)
# ----------------------------
is_poly <- apply(X, 2, function(z) length(unique(z[!is.na(z)])) > 1)
X <- X[, is_poly, drop = FALSE]
if (ncol(X) < 2) stop("Too few polymorphic variants")

miss_rate <- colMeans(is.na(X))
keep_miss <- miss_rate <= miss_cutoff
X <- X[, keep_miss, drop = FALSE]
if (ncol(X) < 2) stop("Too few variants after missing filter")

keep_maf <- NULL
if (maf_cutoff > 0) {
  af <- colMeans(X, na.rm = TRUE)
  maf <- pmin(af, 1 - af)
  keep_maf <- maf >= maf_cutoff
  X <- X[, keep_maf, drop = FALSE]
}

if (impute_missing && anyNA(X)) {
  cm <- colMeans(X, na.rm = TRUE)
  idx <- which(is.na(X), arr.ind = TRUE)
  X[idx] <- cm[idx[,2]]
}

# ----------------------------
# PCA
# ----------------------------
pca <- prcomp(X, center = TRUE, scale. = FALSE)
varprop <- (pca$sdev^2) / sum(pca$sdev^2)
scores <- pca$x

pca_data <- data.frame(
  sample.id = sample_ids,
  PC1 = scores[,1],
  PC2 = scores[,2],
  PC3 = if (ncol(scores) >= 3) scores[,3] else NA,
  PC4 = if (ncol(scores) >= 4) scores[,4] else NA,
  PC5 = if (ncol(scores) >= 5) scores[,5] else NA
)

write.csv(pca_data, file.path(out_dir, "usobtusa_mito_pca.csv"),
          row.names = FALSE, quote = FALSE)

# ----------------------------
# Merge mitotype (NO NA ALLOWED)
# ----------------------------
pca_merged <- merge(
  pca_data, meta,
  by.x = "sample.id",
  by.y = "sample.id",
  all.x = TRUE
)

if (anyNA(pca_merged$mitotype)) stop("ERROR: NA mitotypes in PCA")

write.csv(pca_merged,
          file.path(out_dir, "usobtusa_mito_pca_with_mitotype.csv"),
          row.names = FALSE, quote = FALSE)

# ----------------------------
# Plot
# ----------------------------

subset(pca_merged, is.na(mitotype))
#Gilmer5_H5

p <- ggplot(pca_merged, aes(PC1, PC2, col = mitotype)) +
  geom_point(size = 2) +
  theme_bw() +
  labs(
    title = "mtDNA PCA (filtered samples + variants)",
    x = paste0("PC1 (", round(varprop[1]*100,1), "%)"),
    y = paste0("PC2 (", round(varprop[2]*100,1), "%)")
  )

ggsave(file.path(out_dir, "usobtusa_pca_mito_PC1_PC2_noNA.png"),
       p, width = 7, height = 6, dpi = 300)


ggsave(file.path(out_dir, "usobtusa_pca_mito_PC1_PC2_noNA.pdf"),
       p, width = 7, height = 6, dpi = 300)

# ----------------------------
# Save kept IDs
# ----------------------------
fwrite(data.table(sample.id = sample_ids),
       file.path(out_dir, "kept_samples.txt"),
       col.names = FALSE)

variant_ids <- seqGetData(gds, "variant.id")
variant_ids <- variant_ids[is_poly][keep_miss]
if (maf_cutoff > 0) variant_ids <- variant_ids[keep_maf]

fwrite(data.table(variant.id = variant_ids),
       file.path(out_dir, "kept_variants.txt"),
       col.names = FALSE)

seqClose(gds)
message("DONE")














































































