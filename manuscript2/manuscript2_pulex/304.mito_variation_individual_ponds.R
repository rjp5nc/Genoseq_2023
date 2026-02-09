#ijob -A berglandlab -c10 -p standard --mem=40G
#module load gcc/11.4.0  openmpi/4.1.4 icu R/4.3.1; R

suppressPackageStartupMessages({
  library(SeqArray)
  library(SNPRelate)
  library(ape)
  library(data.table)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(readr)
  library(stringr)
  library(igraph)
  library(ggplot2)
})

# ============================================================
# 0) USER SETTINGS
# ============================================================

LOCATION <- "D8"         # <-- change this to your pond/location label
DIFF_THRESHOLD <- 10        # edges are drawn for diffs < DIFF_THRESHOLD (strictly less than)
MISS_SAMPLE_MAX <- 0.10     # per-sample missingness cutoff (within current sample filter)
MISS_SITE_MAX   <- 0.15     # per-variant missingness cutoff
MIN_MEAN_DP     <- 30       # per-sample mean DP cutoff (computed after filtering)
ROOT_OUTGROUP   <- "RobertUK_A10"  # used if present; otherwise tree is left unrooted

# Base output folder for this run
BASE_OUT <- "/scratch/rjp5nc/UK2022_2024/mitogvcf/gvcf/eudpulex_mito_reverse_out_all"
out_dir  <- file.path(BASE_OUT, LOCATION)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Inputs
gds_path <- "/scratch/rjp5nc/UK2022_2024/redone_mito/eupulex/cohort_gendb/pulex.mito.ALLSITES.rmMissGT0.3.gds"

# Your combined sample table (must include Species, Continent, Sample_ID)
samples_csv <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/Sample_ID_Species_merged_20251227.csv"

# Your SRA merged metadata table
sra_tsv <- "/scratch/rjp5nc/rawdata/sra_metadata_out/sra_merged.tsv"

# Local sequencing metadata with accuratelocation (optional but used for local location filtering)
# If you do not want to use this, set to NA and the script will skip that filter.
metadata_with_clone_csv <- "/project/berglandlab/Robert/UKSequencing2022_2024/2022_2024seqmetadata20250811.csv"

# ============================================================
# 1) HELPERS
# ============================================================

excel_letters <- function(n) {
  # 1 -> A, 26 -> Z, 27 -> AA
  out <- character(n)
  for (i in seq_len(n)) {
    lab <- ""
    x <- i
    while (x > 0) {
      x <- x - 1
      lab <- paste0(LETTERS[(x %% 26) + 1], lab)
      x <- x %/% 26
    }
    out[i] <- lab
  }
  out
}

# Parse SRA SampleName like: March20_2018_DBunk_10  -> date=March20_2018, location=DBunk, clone_id=DBunk_10
parse_sra_samplename <- function(df) {
  df %>%
    tidyr::extract(
      SampleName,
      into = c("date", "location", "clone_num"),
      regex = "^([^_]+_[^_]+)_([^_]+)_(.+)$",
      remove = FALSE
    ) %>%
    dplyr::mutate(
      location = as.character(location),
      clone_id = paste0(location, "_", clone_num)
    )
}

# Given GT matrix (variants x samples) with NA for missing,
# compute pairwise diffs and compared sites
pairwise_diffs <- function(GT) {
  ns <- ncol(GT)
  samp <- colnames(GT)
  is_missing <- is.na(GT) | (GT >= 3L)

  diff_counts <- matrix(0L, ns, ns, dimnames = list(samp, samp))
  compared    <- matrix(0L, ns, ns, dimnames = list(samp, samp))

  for (i in seq_len(ns)) {
    for (j in i:ns) {
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

  list(diff_counts = diff_counts, compared = compared)
}

# Build connected components graph using diffs < threshold
# returns a named vector sample -> group label
assign_groups_from_diffs <- function(diff_counts, compared, threshold) {
  DTd <- as.data.table(as.table(diff_counts))
  setnames(DTd, c("sampleA", "sampleB", "diffs"))
  DTc <- as.data.table(as.table(compared))
  setnames(DTc, c("sampleA", "sampleB", "compared_sites"))
  DT <- merge(DTd, DTc, by = c("sampleA", "sampleB"))
  DT <- DT[sampleA != sampleB]

  edges <- DT[diffs < threshold & compared_sites > 0, .(sampleA, sampleB)]
  g <- igraph::graph_from_data_frame(edges, directed = FALSE, vertices = colnames(diff_counts))
  comp <- igraph::components(g)

  lev <- excel_letters(max(comp$membership))
  group_letters <- setNames(lev[comp$membership], names(comp$membership))
  group_letters
}

# ============================================================
# 2) READ METADATA AND DEFINE THE SAMPLE SET FOR THIS LOCATION
# ============================================================

# Main sample table for species/continent filtering
samples_tbl <- read_csv(samples_csv, show_col_types = FALSE)



# Keep EU D. pulex sample IDs (local names)
eu_pulex_ids <- samples_tbl %>%
  filter(Species == "Daphnia pulex", Continent == "Europe") %>%
  pull(Sample_ID) %>%
  unique()

# Local metadata_with_clone: used to pick samples in this LOCATION (optional)
local_loc_ids <- character(0)
if (!is.na(metadata_with_clone_csv) && file.exists(metadata_with_clone_csv)) {
  meta_loc <- read_csv(metadata_with_clone_csv, show_col_types = FALSE)
  # Common columns you used: clone, Well, accuratelocation
  # We only need Well and accuratelocation
  if (all(c("Well", "accuratelocation") %in% names(meta_loc))) {
    meta_loc <- meta_loc %>%
      mutate(
        Well = trimws(as.character(Well)),
        accuratelocation = trimws(as.character(accuratelocation))
      ) %>%
      filter(!is.na(accuratelocation), accuratelocation != "", accuratelocation == LOCATION) %>%
      distinct(Well) %>%
      pull(Well) %>%
      unique()
    local_loc_ids <- meta_loc
  }
}

# SRA map for D. pulex
sra_map <- read_tsv(sra_tsv, show_col_types = FALSE) %>%
  select(Run, ScientificName, SampleName) %>%
  filter(!is.na(Run), !is.na(SampleName), SampleName != "") %>%
  filter(ScientificName == "Daphnia pulex")

# Parse SRA SampleName into location and clone_id, then keep this LOCATION
sra_parsed <- parse_sra_samplename(sra_map) %>%
  filter(!is.na(location), location == LOCATION)

sra_runs_location <- unique(sra_parsed$Run)

# Combine target samples:
# - EU pulex local IDs
# - local IDs that are in this LOCATION (if available)
# - SRA runs from this LOCATION (these appear in your GDS as SRR...)
target_samples <- unique(c(eu_pulex_ids, local_loc_ids, sra_runs_location))

message("Target samples (union): ", length(target_samples))
message("  local EU pulex IDs: ", length(eu_pulex_ids))
message("  local LOCATION IDs: ", length(local_loc_ids))
message("  SRA runs in LOCATION: ", length(sra_runs_location))

# ============================================================
# 3) OPEN GDS AND FILTER SAMPLES + VARIANTS
# ============================================================

gds <- seqOpen(gds_path)
on.exit(seqClose(gds), add = TRUE)

seqResetFilter(gds)

# Keep only samples present in the GDS
gds_samples <- seqGetData(gds, "sample.id")
target_samples_in_gds <- intersect(target_samples, gds_samples)
if (length(target_samples_in_gds) < 3) stop("Too few target samples found in the GDS for LOCATION = ", LOCATION)

seqSetFilter(gds, sample.id = target_samples_in_gds, verbose = TRUE)

# Per-sample missingness
miss_sample <- seqMissing(gds, per.variant = FALSE)
sample_ids  <- seqGetData(gds, "sample.id")
keep_samp1  <- sample_ids[miss_sample < MISS_SAMPLE_MAX]
if (length(keep_samp1) < 3) stop("Too few samples after sample missingness filter")

seqSetFilter(gds, sample.id = keep_samp1, verbose = TRUE)

# Per-variant missingness filter
miss_var <- seqMissing(gds, per.variant = TRUE)
var_ids  <- seqGetData(gds, "variant.id")
keep_var <- var_ids[miss_var < MISS_SITE_MAX]
if (length(keep_var) < 50) stop("Too few variants after variant missingness filter")

seqSetFilter(gds, sample.id = keep_samp1, variant.id = keep_var, verbose = TRUE)

# Mean DP per sample on the filtered data
DP_raw <- seqGetData(gds, "annotation/format/DP")
# DP_raw commonly variants x samples, but sometimes list
if (is.list(DP_raw)) DP_raw <- do.call(cbind, DP_raw)
DP_mat <- as.matrix(DP_raw)

# Align DP_mat to variants x samples (same as GT after we build it)
# We do not know orientation yet, so compute per sample robustly:
if (ncol(DP_mat) == length(seqGetData(gds, "sample.id"))) {
  mean_dp <- colMeans(DP_mat, na.rm = TRUE)
} else if (nrow(DP_mat) == length(seqGetData(gds, "sample.id"))) {
  mean_dp <- rowMeans(DP_mat, na.rm = TRUE)
} else {
  stop("Cannot infer DP dimensions: ", paste(dim(DP_mat), collapse = " x "))
}

sample_ids2 <- seqGetData(gds, "sample.id")
keep_samp2  <- sample_ids2[!is.na(mean_dp) & mean_dp >= MIN_MEAN_DP]
if (length(keep_samp2) < 3) stop("Too few samples after mean DP filter")

seqSetFilter(gds, sample.id = keep_samp2, variant.id = keep_var, verbose = TRUE)

message("Final kept samples: ", length(seqGetData(gds, "sample.id")))
message("Final kept variants: ", length(seqGetData(gds, "variant.id")))

# ============================================================
# 4) LOAD GENOTYPES (variants x samples)
# ============================================================

GT_raw <- seqGetData(gds, "genotype")  # ploidy x sample x variant OR ploidy x variant x sample
stopifnot(length(dim(GT_raw)) == 3L)

d <- dim(GT_raw)
nsamp <- length(seqGetData(gds, "sample.id"))

# Infer layout
if (d[2] == nsamp) {
  GT <- t(GT_raw[1, , ])        # variants x samples
} else if (d[3] == nsamp) {
  GT <- GT_raw[1, , ]           # variants x samples
} else {
  stop("Cannot infer GT dimension order: dim(genotype) = ", paste(d, collapse = " x "), " ; nsamp = ", nsamp)
}

GT <- as.matrix(GT)
colnames(GT) <- seqGetData(gds, "sample.id")

# Make haploid 0/1 if needed (anything >1 becomes 1)
GT[!is.na(GT) & GT > 1] <- 1L

message("GT matrix: ", nrow(GT), " variants x ", ncol(GT), " samples")

# ============================================================
# 5) PAIRWISE DIFFERENCES (ALL vs ALL) AND DIFFRATE PLOTS
# ============================================================

pw <- pairwise_diffs(GT)
diff_counts <- pw$diff_counts
compared    <- pw$compared
prop_diff   <- diff_counts / compared

# Long table
DT_prop <- as.data.table(as.table(prop_diff)); setnames(DT_prop, c("sampleA","sampleB","prop_diff"))
DT_diff <- as.data.table(as.table(diff_counts)); setnames(DT_diff, c("sampleA","sampleB","diffs"))
DT_comp <- as.data.table(as.table(compared));    setnames(DT_comp, c("sampleA","sampleB","compared_sites"))

pairwise_long <- merge(DT_prop, DT_diff, by = c("sampleA","sampleB"))
pairwise_long <- merge(pairwise_long, DT_comp, by = c("sampleA","sampleB"))
pairwise_long <- pairwise_long[sampleA != sampleB]
pairwise_long[, Similarity := 1 - prop_diff]

fwrite(pairwise_long, file.path(out_dir, "allsite_pairwise_long.tsv"), sep = "\t")

saveRDS(
  list(diff_counts = diff_counts, compared = compared, prop_diff = prop_diff),
  file = file.path(out_dir, "allsite_pairwise_matrices.rds")
)

# Diff-rate histograms
p1 <- ggplot(pairwise_long, aes(x = diffs)) +
  geom_histogram(bins = 100) +
  theme_bw() +
  xlim(0, 100) +
  labs(x = "diffs", y = "Frequency", title = paste0("Pairwise diffs (<=100), ", LOCATION))
ggsave(file.path(out_dir, "diff_rates.png"), p1, width = 7, height = 6, dpi = 300)

p2 <- ggplot(pairwise_long, aes(x = diffs)) +
  geom_histogram(bins = 100) +
  theme_bw() +
  labs(x = "diffs", y = "Frequency", title = paste0("Pairwise diffs (all), ", LOCATION))
ggsave(file.path(out_dir, "diff_rates_all.png"), p2, width = 7, height = 6, dpi = 300)

p3 <- ggplot(pairwise_long, aes(x = diffs)) +
  geom_histogram(bins = 300) +
  theme_bw() +
  xlim(0, 300) +
  labs(x = "diffs", y = "Frequency", title = paste0("Pairwise diffs (<=300), ", LOCATION))
ggsave(file.path(out_dir, "diff_rates_300.png"), p3, width = 7, height = 6, dpi = 300)

# ============================================================
# 6) ASSIGN MITOTYPES FROM DIFF THRESHOLD (< DIFF_THRESHOLD)
# ============================================================

group_letters <- assign_groups_from_diffs(diff_counts, compared, threshold = DIFF_THRESHOLD)

mitotypes_tbl <- tibble(
  sample = names(group_letters),
  mitotype = unname(group_letters)
) %>%
  arrange(mitotype, sample)

write_csv(mitotypes_tbl, file.path(out_dir, "eudpulex_mito_types.csv"))

# Add mitotype into pairwise table (for mitotype-level comparisons)
pairwise_long2 <- pairwise_long %>%
  as_tibble() %>%
  mutate(
    mitoA = mitotypes_tbl$mitotype[match(sampleA, mitotypes_tbl$sample)],
    mitoB = mitotypes_tbl$mitotype[match(sampleB, mitotypes_tbl$sample)]
  )

write_csv(pairwise_long2, file.path(out_dir, "allsite_pairwise_long_with_mitotypes.csv"))

# ============================================================
# 7) PAIRWISE COMPARISONS BETWEEN MITOTYPES (like before)
# ============================================================

# Keep one triangle only to avoid double-counting
pairwise_tri <- pairwise_long2 %>%
  filter(sampleA < sampleB) %>%
  filter(!is.na(mitoA), !is.na(mitoB))

# (A) All pairwise prop_diff grouped by within vs between mitotypes
pairwise_tri <- pairwise_tri %>%
  mutate(comp_type = if_else(mitoA == mitoB, "Within mitotype", "Between mitotypes"))

p_within_between <- ggplot(pairwise_tri, aes(x = comp_type, y = prop_diff)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, height = 0, alpha = 0.4) +
  theme_bw() +
  labs(
    title = paste0("Within vs Between mitotype distances, ", LOCATION),
    x = "",
    y = "Proportion different"
  )
ggsave(file.path(out_dir, "mitotype_within_vs_between_box.png"),
       p_within_between, width = 7, height = 5, dpi = 300)

# (B) Mitotype pair labels (A-B, A-C, etc.), symmetric collapse
pairwise_pairs <- pairwise_tri %>%
  mutate(
    mito_low  = pmin(mitoA, mitoB),
    mito_high = pmax(mitoA, mitoB),
    mito_pair = if_else(mito_low == mito_high,
                        paste0(mito_low, "-", mito_high, " (within)"),
                        paste0(mito_low, "-", mito_high))
  )

p_pairs <- ggplot(pairwise_pairs, aes(x = mito_pair, y = prop_diff)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, height = 0, alpha = 0.35) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1)) +
  labs(
    title = paste0("Pairwise distances by mitotype pair, ", LOCATION),
    x = "Mitotype pair",
    y = "Proportion different"
  )
ggsave(file.path(out_dir, "mitotype_pairwise_box.png"),
       p_pairs, width = 10, height = 6, dpi = 300)

# Optional summary table of mitotype pair means
mitotype_pair_summary <- pairwise_pairs %>%
  group_by(mito_pair) %>%
  summarise(
    n = n(),
    mean_prop_diff = mean(prop_diff, na.rm = TRUE),
    median_prop_diff = median(prop_diff, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(mean_prop_diff)

write_csv(mitotype_pair_summary, file.path(out_dir, "mitotype_pair_summary.csv"))

# ============================================================
# 8) NJ TREE, ROOT, THEN REPLACE SRR TIPS WITH SampleName ONLY
# ============================================================

dist_matrix <- as.dist(prop_diff)
tree <- nj(dist_matrix)

# Root if possible
if (ROOT_OUTGROUP %in% tree$tip.label) {
  tree_rooted <- root(tree, outgroup = ROOT_OUTGROUP, resolve.root = TRUE)
} else {
  tree_rooted <- tree
}

# Clean negative edge lengths
tree_rooted$edge.length[tree_rooted$edge.length < 0] <- 0

# Build SRR -> SampleName map (for this LOCATION only is fine, but safe to keep all)
run_to_sample <- setNames(sra_map$SampleName, sra_map$Run)

# Drop SRR tips with no mapping so final plot has NO SRR labels
unmapped_srr <- tree_rooted$tip.label[
  grepl("^SRR", tree_rooted$tip.label) & !(tree_rooted$tip.label %in% names(run_to_sample))
]
if (length(unmapped_srr) > 0) {
  tree_rooted <- drop.tip(tree_rooted, unmapped_srr)
}

# Replace SRR -> SampleName
tip_new <- tree_rooted$tip.label
is_srr  <- tip_new %in% names(run_to_sample)
tip_new[is_srr] <- unname(run_to_sample[tip_new[is_srr]])

# Keep everything else unchanged
tree_rooted$tip.label <- tip_new

# Ensure unique tips
tree_rooted$tip.label <- make.unique(tree_rooted$tip.label, sep = "_dup")

# Save tree
write.tree(tree_rooted, file = file.path(out_dir, "eudpulex_tree_nwk.nwk"))

# ============================================================
# 9) COLOR TREE BY MITOTYPE (mapped into SampleName space)
# ============================================================

# Convert mitotypes table into SampleName space:
# If a "sample" is SRR and exists in run_to_sample, map it to SampleName.
mitotypes_tbl2 <- mitotypes_tbl %>%
  mutate(
    sample2 = if_else(sample %in% names(run_to_sample),
                      unname(run_to_sample[sample]),
                      sample)
  ) %>%
  select(sample = sample2, mitotype) %>%
  mutate(sample = make.unique(sample, sep = "_dup"))

mitotype_map2 <- setNames(mitotypes_tbl2$mitotype, mitotypes_tbl2$sample)
tip_mito <- unname(mitotype_map2[tree_rooted$tip.label])

mito_levels <- sort(unique(na.omit(tip_mito)))
mito_cols   <- setNames(rainbow(length(mito_levels)), mito_levels)

tip_colors <- mito_cols[tip_mito]
tip_colors[is.na(tip_colors)] <- "grey50"

xmax <- as.numeric(quantile(node.depth.edgelength(tree_rooted), 0.98))

# Circular PDF
pdf(file.path(out_dir, "eudpulex_groups_circ.pdf"), width = 20, height = 15)
plot.phylo(
  tree_rooted,
  type = "fan",
  cex = 0.7,
  no.margin = TRUE,
  tip.color = tip_colors,
  x.lim = c(0, xmax),
  main = paste0("Rooted NJ Tree (mito, all sites), ", LOCATION)
)
legend(
  "topleft",
  legend = mito_levels,
  col = mito_cols[mito_levels],
  pch = 19,
  pt.cex = 1.3,
  cex = 0.9,
  bty = "n",
  title = "Mitotype"
)
dev.off()

# Circular PNG
png(file.path(out_dir, "eudpulex_groups_circ.png"), width = 4500, height = 3500, res = 300)
plot.phylo(
  tree_rooted,
  type = "fan",
  cex = 0.7,
  no.margin = TRUE,
  tip.color = tip_colors,
  x.lim = c(0, xmax),
  main = paste0("Rooted NJ Tree (mito, all sites), ", LOCATION)
)
legend(
  "topleft",
  legend = mito_levels,
  col = mito_cols[mito_levels],
  pch = 19,
  pt.cex = 1.3,
  cex = 0.9,
  bty = "n",
  title = "Mitotype"
)
dev.off()

# ============================================================
# 10) BOXPLOTS WITHIN THIS LOCATION (mitotypes within pond)
# ============================================================

# If you truly only want within-LOCATION comparisons, restrict pairwise rows to samples whose
# displayed names indicate LOCATION for SRA and local IDs.
# For SRR tips, we already mapped to SampleName, so filter by SampleName containing "_<LOCATION>_"
# Example SampleName: March20_2018_DBunk_23  contains "_DBunk_"
loc_pat <- paste0("_", LOCATION, "_")

pairwise_loc <- pairwise_long2 %>%
  mutate(
    sampleA2 = if_else(sampleA %in% names(run_to_sample), unname(run_to_sample[sampleA]), sampleA),
    sampleB2 = if_else(sampleB %in% names(run_to_sample), unname(run_to_sample[sampleB]), sampleB)
  ) %>%
  filter(str_detect(sampleA2, fixed(loc_pat)) | sampleA2 %in% local_loc_ids | sampleA2 %in% eu_pulex_ids) %>%
  filter(str_detect(sampleB2, fixed(loc_pat)) | sampleB2 %in% local_loc_ids | sampleB2 %in% eu_pulex_ids) %>%
  mutate(
    mitoA = mitotypes_tbl2$mitotype[match(make.unique(sampleA2, sep = "_dup"), mitotypes_tbl2$sample)],
    mitoB = mitotypes_tbl2$mitotype[match(make.unique(sampleB2, sep = "_dup"), mitotypes_tbl2$sample)]
  )

pairwise_loc_tri <- pairwise_loc %>%
  filter(sampleA < sampleB) %>%
  filter(!is.na(mitoA), !is.na(mitoB)) %>%
  mutate(
    mito_low  = pmin(mitoA, mitoB),
    mito_high = pmax(mitoA, mitoB),
    mito_pair = if_else(mito_low == mito_high,
                        paste0(mito_low, "-", mito_high, " (within)"),
                        paste0(mito_low, "-", mito_high))
  )

p_loc_pairs <- ggplot(pairwise_loc_tri, aes(x = mito_pair, y = prop_diff)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, height = 0, alpha = 0.35) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1)) +
  labs(
    title = paste0("Within-", LOCATION, " distances by mitotype pair"),
    x = "Mitotype pair",
    y = "Proportion different"
  )

ggsave(file.path(out_dir, paste0("within_", LOCATION, "_mitotype_pairwise_box.png")),
       p_loc_pairs, width = 10, height = 6, dpi = 300)

# ============================================================
# 11) SAVE KEPT SAMPLE LIST (as present in GDS)
# ============================================================

writeLines(seqGetData(gds, "sample.id"), file.path(out_dir, "kept_samples.txt"))
writeLines(as.character(seqGetData(gds, "variant.id")), file.path(out_dir, "kept_variants.txt"))

message("DONE. Outputs in: ", out_dir)


# DP on the CURRENT filtered GDS (same samples/variants used above)
DP_raw2 <- seqGetData(gds, "annotation/format/DP")
if (is.list(DP_raw2)) DP_raw2 <- do.call(cbind, DP_raw2)
DP_mat2 <- as.matrix(DP_raw2)

kept_ids_now <- seqGetData(gds, "sample.id")

# Compute mean DP per sample robustly regardless of DP orientation
mean_dp_now <- if (ncol(DP_mat2) == length(kept_ids_now)) {
  colMeans(DP_mat2, na.rm = TRUE)
} else if (nrow(DP_mat2) == length(kept_ids_now)) {
  rowMeans(DP_mat2, na.rm = TRUE)
} else {
  stop("Cannot infer DP matrix orientation: dim(DP) = ",
       paste(dim(DP_mat2), collapse = " x "),
       " ; n_samples = ", length(kept_ids_now))
}

dp_tbl <- tibble(
  sample = kept_ids_now,
  meanDP = as.numeric(mean_dp_now)
)

# Attach mitotype (SRR IDs are stored as SRR in mitotypes_tbl;
# local IDs remain unchanged)
dp_mito <- dp_tbl %>%
  left_join(mitotypes_tbl, by = "sample") %>%
  filter(!is.na(mitotype))

if (nrow(dp_mito) == 0) stop("No samples with mitotype assignment found when selecting representatives.")

# One sample per mitotype: highest DP (tie-breaker = alphabetical sample)
rep_per_mitotype <- dp_mito %>%
  arrange(mitotype, desc(meanDP), sample) %>%
  group_by(mitotype) %>%
  slice(1) %>%
  ungroup()

# Save a table + a plain sample list
write_csv(rep_per_mitotype, file.path(out_dir, "representative_per_mitotype_highestDP.csv"))
writeLines(rep_per_mitotype$sample, file.path(out_dir, "representative_per_mitotype_highestDP.txt"))

# Optional: if you want the representative names in SampleName space (SRR -> SampleName),
# save an extra file with "display_name"
rep_per_mitotype_display <- rep_per_mitotype %>%
  mutate(
    display_name = if_else(sample %in% names(run_to_sample),
                           unname(run_to_sample[sample]),
                           sample)
  )

write_csv(rep_per_mitotype_display,
          file.path(out_dir, "representative_per_mitotype_highestDP_withSampleName.csv"))

message("Saved representatives: ",
        file.path(out_dir, "representative_per_mitotype_highestDP.csv"))