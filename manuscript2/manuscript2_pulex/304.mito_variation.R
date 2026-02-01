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

library(patchwork)
library(foreach)
library(lubridate)




# ---- Step 1: Open GDS file ----

metadata <- read.csv("/project/berglandlab/Robert/UKSequencing2022_2024/old_stuff/2022_2024seqmetadata20250131.csv", header = TRUE)
metadata_with_clone <- read.csv("/project/berglandlab/Robert/UKSequencing2022_2024/old_stuff/metadata_with_clone.csv", header = TRUE)

samples <- read.csv("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/Sample_ID_Species_merged_20251227.csv")

usobtusasamps <- subset(samples, Species == "Daphnia obtusa" & Continent == "NorthAmerica")

metadata_with_clone <- read.csv("/project/berglandlab/Robert/UKSequencing2022_2024/old_stuff/metadata_with_clone.csv", header = TRUE)
metadata_with_clone$clone <- trimws(metadata_with_clone$clone)
metadata_with_clone <- subset(metadata_with_clone, clone !="Blank")
metadata_with_clone <- subset(metadata_with_clone, clone !="BLANK")

diffs <- read.csv("/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/mito_types_diff_df.csv")


head(metadata)
metadata$clone <- trimws(metadata$clone)
metadata <- metadata %>% 
  filter(!tolower(clone) %in% c("blank", "blanks", "na", "missing"))

seqClose(genofile)
#gds.fn <- "/scratch/rjp5nc/UK2022_2024/mito_vcf/usdobtusa_mito_genotyped.gds"             # output GDS file
gds.fn <- "/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/usdobtusa_mito_allsites_all.haploid.annotated2.gds"
genofile <- seqOpen(gds.fn)

seqResetFilter(genofile)

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

seqSetFilter(genofile, sample.id = final_valid_samples)

miss_rate <- seqMissing(genofile, per.variant = TRUE)
dp <- seqGetData(genofile, "annotation/format/DP")
mean_depth <- rowMeans(dp, na.rm = TRUE)
keep <- which(miss_rate < 0.10)

g_samples <- seqGetData(genofile, "sample.id")

mask <- (mean_depth > 40) & !is.na(mean_depth)
keep_samples <- g_samples[mask]
keep_samples <- keep_samples[!is.na(keep_samples)]  # drop any NA slots
stopifnot(all(!is.na(keep_samples)))
final_valid_samples <- intersect(keep_samples, final_valid_samples)

samps300 <- subset(diffs, diffs < 500)
samps300samps <- samps300$sample

final_valid_samples <- intersect(samps300samps, final_valid_samples)


# ---- Step 3: Compute IBS distance matrix ----
ibs <- snpgdsIBS(genofile, sample.id = final_valid_samples, snp.id = keep, num.thread = 4, autosome.only = FALSE)
dist_matrix <- 1 - ibs$ibs
rownames(dist_matrix) <- colnames(dist_matrix) <- ibs$sample.id
write.csv(dist_matrix, "/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/usobtusa_mito_similarity_matrix.csv")  # Save as CSV

dist_matrix2 <- 1- dist_matrix

long_mat <- melt(dist_matrix2, varnames = c("CloneA", "CloneB"), value.name = "Similarity") %>%
  mutate(CloneA = as.character(CloneA),
         CloneB = as.character(CloneB))  # ensure they are character, not factor

long_mat



long_mat_plot <- ggplot(long_mat, aes(x = Similarity)) +
  geom_histogram(bins = 200) +
  theme_bw() + 
  xlab("Similarity") + ylim(0,1000000) + xlim(0.5,1.05) + 
  ylab("Frequency") +
  geom_vline(xintercept = 0.975, linetype = "dashed", color = "red") 

ggsave("/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/long_mat_mito_plot.png", plot = long_mat_plot, width = 6, height = 4, dpi = 300)


long_mat_plot <- ggplot(long_mat, aes(x = Similarity)) +
  geom_histogram(bins = 200) +
  theme_bw() + 
  xlab("Similarity") +
  xlim(0.5, 1.05) +
  ylab("Frequency") +
  scale_y_log10() +
  geom_vline(xintercept = 0.975, linetype = "dashed", color = "red")
ggsave("/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/long_mat_mito_plot.png", plot = long_mat_plot, width = 6, height = 4, dpi = 300)


options(digits = 16)
head(long_mat$Similarity)

# 2. Filter dataset
long_filt <- long_mat %>%
  filter(Similarity >= 0.9999 & CloneA != CloneB)

#

# 3. Build graph
g <- graph_from_data_frame(long_filt[, c("CloneA", "CloneB")], directed = FALSE)

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
long_filt <- long_filt %>%
  mutate(Group = group_letters[CloneA])


head(long_filt)
unique(long_filt$Group)

long_filt_onlywell <- long_filt[,c(1,4)]

unique_clones <- long_filt_onlywell %>%
  distinct(CloneA, Group)
unique_clones
unique(long_filt$Group)


# ---- Step 4: Build neighbor-joining tree ----
tree <- nj(as.dist(dist_matrix))

# ---- Step 5 (Optional): Relabel tips from metadata CSV ----
# CSV should have columns: Sample_ID and Label (or Well, Clone, etc.)
tree_rooted <- root(tree, outgroup = "Rockpool2_C6", resolve.root = TRUE)

label_map <- setNames(unique_clones$CloneA, unique_clones$Group)

# Apply labels (optional: fallback to Sample_ID if label is missing)
tree$tip.label <- label_map[tree$tip.label]
tree$tip.label[is.na(tree$tip.label)] <- names(label_map)[is.na(tree$tip.label)]

group_colors <- setNames(rainbow(length(unique(unique_clones$Group))), unique(unique_clones$Group))

# Assign color to each tip

tip_colors <- group_colors[unique_clones$Group]

tip_color_dt <- data.table(unique_clones$Group, unique_clones$CloneA)

final_valid_samples2 <- as.data.frame(final_valid_samples)

final_valid_samples3 <- left_join(final_valid_samples2, tip_color_dt, by= c("final_valid_samples"="V2"))

tip_colors <- group_colors[final_valid_samples3$V1]


tip_color_dt <- data.table(
  sample = names(tip_colors),
  color  = as.vector(tip_colors)
)


png("/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/usobtusa_mito_groups.png", res = 300, width = 4000, height = 9000)
#png("/scratch/rjp5nc/UK2022_2024/mito_vcf/tree_usobtusa_circ.png", width = 1200, height = 2000)

# Plot the tree
plot.phylo(tree_rooted,
           type = "phylogram",
#           type = "fan",
           cex = 0.8,
           label.offset = 0.01,
           no.margin = TRUE,
           tip.color = tip_colors,
           main = "Rooted NJ Tree (root = P759)")

# Add legend
legend("topleft",                   # position
       legend = unique(tip_color_dt$sample),     # group names
       col = unique(tip_color_dt$color),        # matching colors
       pch = 19,                   # solid circle
       pt.cex = 1.5,               # point size
       cex = 1,                    # text size
       bty = "n",                  # no box
       title = "Sample Group")

dev.off()

write.tree(tree_rooted, file = "/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/usdobtusa_tree_nwk.nwk")

seqClose(genofile)

write.csv(unique_clones, "/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/mito_types.csv")  # Save as CSV




gds_path <- "/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/usdobtusa_mito_allsites_all.haploid.annotated2.gds"
mitotype_csv <- "/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/mito_types.csv"
out_dir <- "/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/unique_snps_by_group"
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

mask <- (mean_dp > 40) & !is.na(mean_dp)
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
mitotypes <- read.csv("/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/mito_types_v2.csv",
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

mitotypes <- subset(mitotypes, Group %in% c("A", "B", "C", "D", "E", "F", "XX"))
mitotypesE <- subset(mitotypes, Group %in% c("E", "F", "XX"))

mito_keep <- mitotypes %>%
  transmute(sample = trimws(as.character(CloneA)),
            group  = trimws(as.character(Group))) %>%
  semi_join(tibble(sample = samps), by = "sample")


write.table(mito_keep$sample, "/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/alignment_us/keep_samplesdp20.txt",
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)



write.table(mito_keep, "/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/alignment_us/conversion.txt",
            sep = "\t",
            row.names = TRUE,
            quote = FALSE)


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


ref_sample <-  "Gilmer5_B4"
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
out_dir <- "/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/unique_snps_by_group"
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
       p_box, width = 6, height = 5, dpi = 600)

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


diff_df


hist <- ggplot(diff_df, aes(x = diffs)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "black") +
  labs(
    title = "Distribution of callable sites per window",
    x = "Callable sites per window",
    y = "Count of windows"
  ) +
  theme_bw(base_size = 12)

ggsave(file.path(out_dir, "hist.png"),
       hist, width = 10, height = 4, dpi = 300)

















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
out_dir <- "/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/unique_snps_by_group"
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
























bed <- read.delim(
  "/scratch/rjp5nc/Reference_genomes/mito_reference/usdobtusa_annotation/usdobtusa.bed",
  header = FALSE, sep = "\t",
  col.names = c("chr","start","end","name","score","strand"),
  stringsAsFactors = FALSE
)

# Feature type (for optional styling/labeling)
bed$type <- ifelse(grepl("^trn", bed$name), "tRNA",
                   ifelse(grepl("^rrn", bed$name), "rRNA", "protein_coding"))
# Only label non-tRNA features to reduce overlap
bed$label <- ifelse(bed$type == "tRNA", NA, bed$name)

# Match your facet order
bed$chr <- factor(bed$chr, levels = chr_levels)

# ---- Helper to add annotations to an existing plot ----
add_mito_annots <- function(p, ann = bed) {
  p +
    # Light bands spanning the y-range for each annotated interval
    geom_rect(
      data = ann,
      aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE,
      fill = "grey60", alpha = 0.15, show.legend = FALSE
    ) +
    # Gene labels at the top margin (non-tRNA only)
    geom_text(
      data = subset(ann, !is.na(label)),
      aes(x = (start + end) / 2, y = Inf, label = label),
      inherit.aes = FALSE, vjust = 1.1, size = 2.6
    ) +
    # Let labels draw outside the panel a bit
    coord_cartesian(clip = "off") +
    theme(plot.margin = margin(14, 12, 14, 12))
}

# ---- Apply to your plots and save ----
p_tracks_ann <- add_mito_annots(p_tracks)
p_group_ann  <- add_mito_annots(p_group)

ggsave(file.path(out_dir, paste0("tracks_across_genome_vs_",
       make.names(ref_sample), "_", window_size/1000, "kb.with_annots.png")),
       p_tracks_ann, width = 10, height = max(4, length(chr_levels) * 1.2), dpi = 300)

ggsave(file.path(out_dir, paste0("groupmean_across_genome_vs_",
       make.names(ref_sample), "_", window_size/1000, "kb.with_annots.png")),
       p_group_ann, width = 10, height = max(4, length(chr_levels) * 1.2), dpi = 300)




























suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(GenomicRanges)
})

# --- BED (already loaded earlier) ---
# Ensure genomic order
bed <- bed %>%
  arrange(chr, start) %>%
  mutate(
    type = case_when(
      grepl("^trn", name) ~ "tRNA",
      grepl("^rrn", name) ~ "rRNA",
      TRUE ~ "protein_coding"
    ),
    gene = name
  )

# Factor level order = genomic order
gene_levels <- bed$gene

# --- Build GRanges for BED ---
gr_bed <- GRanges(
  seqnames = bed$chr,
  ranges   = IRanges(start = bed$start + 1, end = bed$end),
  gene     = bed$gene,
  type     = bed$type
)

window_size = 1

# --- Windows (use existing start/end if present; else derive from mid & window_size) ---
if (!all(c("start","end") %in% names(win_df))) {
  win_df <- win_df %>%
    mutate(
      start = pmax(1L, floor(mid - window_size / 2)),
      end   = ceiling(mid + window_size / 2)
    )
}

# Keep chr order consistent
win_df$chr <- factor(win_df$chr, levels = chr_levels)

gr_win <- GRanges(
  seqnames = win_df$chr,
  ranges   = IRanges(start = win_df$start, end = win_df$end)
)

# --- Assign each window to the gene with MAX overlap ---
hits <- findOverlaps(gr_win, gr_bed)
if (length(hits) == 0L) stop("No window ↔ gene overlaps found. Check chr names and coordinates.")

ov <- tibble(
  wid  = queryHits(hits),
  gid  = subjectHits(hits),
  ovw  = width(pintersect(gr_win[queryHits(hits)], gr_bed[subjectHits(hits)]))
) %>%
  group_by(wid) %>% slice_max(ovw, n = 1, with_ties = FALSE) %>% ungroup()

joined <- bind_cols(
  win_df[ov$wid, , drop = FALSE],
  as.data.frame(mcols(gr_bed[ov$gid]))
)

# Order gene factor by genome order
joined$gene <- factor(joined$gene, levels = gene_levels)

# --- Data frames for plotting ---
df_all    <- joined
df_notrna <- joined %>% filter(type != "tRNA")

# ------------------ PLOTS ------------------

# A) All features
p_boxes_all <- ggplot(df_all, aes(x = gene, y = prop_diff)) +
  geom_boxplot(outlier.size = 0.4) +
  facet_wrap(~ chr, scales = "free_x", ncol = 1, strip.position = "left") +
  labs(
    title = paste0("Differences vs ", ref_sample, " by gene region (", window_size/1000, " kb windows)"),
    x = "Gene / Feature",
    y = "Proportion different"
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.spacing.y   = unit(0.4, "lines"),
    strip.background  = element_rect(fill = "grey90"),
    strip.placement   = "outside",
    axis.text.x       = element_text(angle = 45, hjust = 1)
  )

ggsave(file.path(out_dir, paste0("boxplots_by_gene_all_vs_",
       make.names(ref_sample), "_", window_size/1000, "kb.png")),
       p_boxes_all, width = 12, height = max(4, length(chr_levels) * 2), dpi = 300)

# B) Non-tRNA
p_boxes_notrna <- ggplot(df_notrna, aes(x = gene, y = prop_diff)) +
  geom_boxplot(outlier.size = 0.4) +
  facet_wrap(~ chr, scales = "free_x", ncol = 1, strip.position = "left") +
  labs(
    title = paste0("Differences vs ", ref_sample, " by gene region (non-tRNA; ",
                   window_size/1000, " kb windows)"),
    x = "Gene",
    y = "Proportion different"
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.spacing.y   = unit(0.4, "lines"),
    strip.background  = element_rect(fill = "grey90"),
    strip.placement   = "outside",
    axis.text.x       = element_text(angle = 45, hjust = 1)
  )

ggsave(file.path(out_dir, paste0("boxplots_by_gene_notRNA_vs_",
       make.names(ref_sample), "_", window_size/1000, "kb.png")),
       p_boxes_notrna, width = 12, height = max(4, length(chr_levels) * 2), dpi = 300)

# C) Optional: by Group (only if present)
if ("group" %in% names(df_notrna)) {
  p_boxes_bygroup <- ggplot(df_notrna, aes(x = gene, y = prop_diff, fill = group)) +
    geom_boxplot(outlier.size = 0.3, position = position_dodge(width = 0.75)) +
    facet_wrap(group ~ chr, scales = "free_x", ncol = 1, strip.position = "left") +
    labs(
      title = paste0("Differences vs ", ref_sample, " by gene region and group (non-tRNA)"),
      x = "Gene",
      y = "Proportion different",
      fill = "Group"
    ) +
    theme_bw(base_size = 10) +
    theme(
      panel.spacing.y   = unit(0.4, "lines"),
      strip.background  = element_rect(fill = "grey90"),
      strip.placement   = "outside",
      axis.text.x       = element_text(angle = 45, hjust = 1)
    )

  ggsave(file.path(out_dir, paste0("boxplots_by_gene_notRNA_bygroup_vs_",
         make.names(ref_sample), "_", window_size/1000, "kb.png")),
         p_boxes_bygroup, width = 12, height = max(4, length(chr_levels) * 15), dpi = 300)
}




















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















#PCA




suppressPackageStartupMessages({
  library(SeqArray)
  library(data.table)
  library(ggplot2)
})

# ----------------------------
# Paths
# ----------------------------
mitotype_csv <- "/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/mito_types_v2.csv"
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
setnames(meta, c("idx", "sample.id", "mitotype"))

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
keep_samples <- intersect(gds_samples, meta_samples)
seqSetFilter(gds, sample.id = keep_samples, verbose = TRUE)

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
keep_samples <- intersect(keep_samples, meta_samples)
seqSetFilter(gds, sample.id = keep_samples, verbose = TRUE)

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
pca_merged <- merge(pca_data, meta, by = "sample.id", all.x = TRUE)
if (anyNA(pca_merged$mitotype)) stop("ERROR: NA mitotypes in PCA")

write.csv(pca_merged,
          file.path(out_dir, "usobtusa_mito_pca_with_mitotype.csv"),
          row.names = FALSE, quote = FALSE)

# ----------------------------
# Plot
# ----------------------------
p <- ggplot(pca_merged, aes(PC1, PC2, col = mitotype)) +
  geom_point(size = 2) +
  theme_bw() +
  labs(
    title = "mtDNA PCA (filtered samples + variants)",
    x = paste0("PC1 (", round(varprop[1]*100,1), "%)"),
    y = paste0("PC2 (", round(varprop[2]*100,1), "%)")
  )

ggsave(file.path(out_dir, "usobtusa_pca_mito_PC1_PC2.png"),
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



# PCA














#Keep NA's





suppressPackageStartupMessages({
  library(SeqArray)
  library(data.table)
  library(ggplot2)
})

# ----------------------------
# Paths
# ----------------------------
mitotype_csv <- "/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/mito_types_v2.csv"
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
setnames(meta, c("idx", "sample.id", "mitotype"))

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
seqSetFilter(gds, sample.id = keep_samples, verbose = TRUE)

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
pca_merged <- merge(pca_data, meta, by = "sample.id", all.x = TRUE)
if (anyNA(pca_merged$mitotype)) stop("ERROR: NA mitotypes in PCA")

write.csv(pca_merged,
          file.path(out_dir, "usobtusa_mito_pca_with_mitotype.csv"),
          row.names = FALSE, quote = FALSE)

# ----------------------------
# Plot
# ----------------------------
p <- ggplot(pca_merged, aes(PC1, PC2, col = mitotype)) +
  geom_point(size = 2) +
  theme_bw() +
  labs(
    title = "mtDNA PCA (filtered samples + variants)",
    x = paste0("PC1 (", round(varprop[1]*100,1), "%)"),
    y = paste0("PC2 (", round(varprop[2]*100,1), "%)")
  )

ggsave(file.path(out_dir, "usobtusa_pca_mito_PC1_PC2.png"),
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



# PCA




















suppressPackageStartupMessages({
  library(SeqArray)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(data.table)
  library(readr)
})

# ----------------------------
# Assumes you already opened/filtered the GDS and have:
#   g         = SeqArray/SeqVarGDSClass object (already subset to samples/variants you want)
#   samps     = character vector of sample IDs in the current GDS selection (seqGetData(g,"sample.id"))
#   mito_keep = data.frame/tibble with at least: sample, group  (optional but used for joins)
# ----------------------------

# ----------------------------
# 4) Pull ALL sites genotypes; ensure GT is [variants x samples]
# ----------------------------

seqResetFilter(genofile)
seqSetFilter(genofile, sample.id = final_valid_samples, variant.id = keep)


GT_raw <- seqGetData(genofile, "genotype")  # haploid: ploidy x variants x samples
if (length(dim(GT_raw)) == 3L) {
  GT <- matrix(GT_raw[1, , ], nrow = dim(GT_raw)[2], ncol = dim(GT_raw)[3])
} else if (is.matrix(GT_raw)) {
  GT <- GT_raw
} else stop("Unexpected GT shape.")

# Fix orientation if needed
if (ncol(GT) != length(samps) && nrow(GT) == length(samps)) GT <- t(GT)
stopifnot(ncol(GT) == length(samps))

# Missing calls: NA or codes >= 3 (adjust if your encoding differs)
is_missing <- is.na(GT) | (GT >= 3L)
ok_call <- !is_missing

# ----------------------------
# 6) All-vs-all pairwise differences across ALL sites
#    Outputs:
#      diff_counts[i,j] = #sites both called and different
#      compared[i,j]    = #sites both called
#      prop_diff[i,j]   = diff_counts/compared
# ----------------------------
nsamp <- length(samps)

diff_counts <- matrix(0L, nsamp, nsamp, dimnames = list(samps, samps))
compared    <- matrix(0L, nsamp, nsamp, dimnames = list(samps, samps))

for (i in seq_len(nsamp)) {
  gi <- GT[, i]
  oi <- ok_call[, i]

  for (j in i:nsamp) {
    ok <- oi & ok_call[, j]
    n_ok <- sum(ok)

    compared[i, j] <- n_ok
    compared[j, i] <- n_ok

    if (n_ok > 0) {
      d <- sum(gi[ok] != GT[ok, j])
      diff_counts[i, j] <- d
      diff_counts[j, i] <- d
    }
  }
}

prop_diff <- diff_counts / compared
prop_diff[!is.finite(prop_diff)] <- NA_real_

diag(diff_counts) <- 0L
diag(prop_diff) <- 0
# Optional: compared on diagonal as per-sample called sites
diag(compared) <- as.integer(colSums(ok_call))

# ----------------------------
# 7) Long table for plotting / thresholding
# ----------------------------
pairwise_long <- as.data.table(as.table(prop_diff))
setnames(pairwise_long, c("sampleA", "sampleB", "prop_diff"))

pairwise_long[, diffs := as.vector(diff_counts)]
pairwise_long[, compared_sites := as.vector(compared)]
pairwise_long[, Similarity := 1 - prop_diff]

# Remove self-comparisons
pairwise_long <- pairwise_long[sampleA != sampleB]

# Attach group labels (if mito_keep present)
if (exists("mito_keep") && all(c("sample", "group") %in% colnames(mito_keep))) {
  pairwise_long <- pairwise_long %>%
    left_join(mito_keep %>% select(sample, group), by = c("sampleA" = "sample")) %>%
    rename(groupA = group) %>%
    left_join(mito_keep %>% select(sample, group), by = c("sampleB" = "sample")) %>%
    rename(groupB = group)
}

# ----------------------------
# 8) Save outputs
# ----------------------------
out_dir <- "/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/unique_snps_by_group"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

saveRDS(
  list(diff_counts = diff_counts, compared = compared, prop_diff = prop_diff),
  file = file.path(out_dir, "allsite_pairwise_matrices.rds")
)

write_tsv(pairwise_long, file.path(out_dir, "allsite_pairwise_long.tsv"))

# Also write a compact CSV if you want
fwrite(pairwise_long, file.path(out_dir, "allsite_pairwise_long.csv"))

message("Wrote:\n  ", file.path(out_dir, "allsite_pairwise_matrices.rds"),
        "\n  ", file.path(out_dir, "allsite_pairwise_long.tsv"),
        "\n  ", file.path(out_dir, "allsite_pairwise_long.csv"))

# ----------------------------
# 9) Optional: quick histogram of similarity (log y)
# ----------------------------
p_hist <- ggplot(pairwise_long, aes(x = Similarity)) +
  geom_histogram(bins = 200) +
  theme_bw() +
  xlab("Similarity (1 - prop_diff)") +
  ylab("Frequency") +
  scale_y_log10() +
  geom_vline(xintercept = 0.9995, linetype = "dashed", color = "red")

ggsave(file.path(out_dir, "allsite_similarity_hist_logy.png"),
       p_hist, width = 7, height = 5, dpi = 300)

# Return objects for interactive use
list(
  diff_counts = diff_counts,
  compared = compared,
  prop_diff = prop_diff,
  pairwise_long = pairwise_long
)


summary(pairwise_long$compared_sites)


