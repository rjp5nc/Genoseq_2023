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
# 0) GLOBAL SETTINGS
# ============================================================

DIFF_THRESHOLD  <- 7
MISS_SAMPLE_MAX <- 0.10
MISS_SITE_MAX   <- 0.15
MIN_MEAN_DP     <- 30
ROOT_OUTGROUP   <- "RobertUK_H6"

BASE_OUT <- "/scratch/rjp5nc/UK2022_2024/mitogvcf/gvcf/eudobtusa_mito_reverse_out_all"
gds_path <- "/scratch/rjp5nc/UK2022_2024/redone_mito/euobtusa/cohort_gendb/obtusa.mito.ALLSITES.rmMissGT0.3.gds"
samples_csv <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/Sample_ID_Species_merged_20251227.csv"
sra_tsv <- "/scratch/rjp5nc/rawdata/sra_metadata_out/sra_merged.tsv"
metadata_with_clone_csv <- "/project/berglandlab/Robert/UKSequencing2022_2024/2022_2024seqmetadata20250811.csv"

# >>> SET THIS ONCE <<<
# Specify ONE or MULTIPLE pond locations (case-insensitive, spaces ok)


LOCATIONS <- c(
  "beavercreek", "birdhut", "cafe", "canal",
  "diris", "dlily", "dmountie", "dmud", "dnorden", "doak", "doily",
  "drail", "dramps", "drusty", "dsandy", "dtiny", "elvis", "islands",
  "melford", "p17", "p19", "p25", "p63", "pbo66", "p30", "p6",
  "pottcarr1", "pottcarr2", "pottcarr3", "psnake", "pvet",
  "shallow", "sharp", "warbler", "xxx5"
)


# If you want: LOCATIONS_STR <- "Dcat, dbUnk, FS"; LOCATIONS <- strsplit(LOCATIONS_STR, ",")[[1]]

# Overlap control (prevents empty tree / empty mitotypes)
MIN_PAIRWISE_COMPARED <- 200   # require >= this many overlapping sites per pair to consider "connected"
MIN_OK_PAIRS_PER_SAMPLE <- 2   # each kept sample must have >= this many partners with overlap >= MIN_PAIRWISE_COMPARED
MIN_SAMPLES_FOR_TREE <- 3      # need at least 3

# ============================================================
# 1) HELPERS
# ============================================================

excel_letters <- function(n) {
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

# IMPORTANT: handle SeqArray missing genotype encodings
pairwise_diffs <- function(GT) {
  ns <- ncol(GT)
  samp <- colnames(GT)

  # Missing often appears as -1 in SeqArray genotype arrays
  is_missing <- is.na(GT) | (GT < 0L) | (GT >= 3L)

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
  setNames(lev[comp$membership], names(comp$membership))
}

mean_dp_per_sample <- function(DP_raw, sample_ids) {
  if (is.list(DP_raw)) DP_raw <- do.call(cbind, DP_raw)
  DP_mat <- as.matrix(DP_raw)

  if (ncol(DP_mat) == length(sample_ids)) {
    return(colMeans(DP_mat, na.rm = TRUE))
  } else if (nrow(DP_mat) == length(sample_ids)) {
    return(rowMeans(DP_mat, na.rm = TRUE))
  } else {
    stop("Cannot infer DP matrix orientation: dim(DP)=",
         paste(dim(DP_mat), collapse = "x"),
         " ; n_samples=", length(sample_ids))
  }
}

canon_loc <- function(x) {
  x <- trimws(as.character(x))
  tolower(x)
}

clone_prefix_before_underscore <- function(x) {
  x <- trimws(as.character(x))
  has_us <- grepl("_", x)
  out <- rep(NA_character_, length(x))
  out[has_us] <- sub("_.*$", "", x[has_us])
  out
}

parse_sra_location_from_samplename <- function(sample_name) {
  m <- str_match(sample_name, "^([^_]+_[^_]+)_([^_]+)_.+$")
  m[, 3]
}

make_tree_tiplabels_sample_only <- function(tree, run_to_sample) {
  unmapped_srr <- tree$tip.label[
    grepl("^SRR", tree$tip.label) & !(tree$tip.label %in% names(run_to_sample))
  ]
  if (length(unmapped_srr) > 0) tree <- drop.tip(tree, unmapped_srr)

  tip_new <- tree$tip.label
  is_srr <- tip_new %in% names(run_to_sample)
  tip_new[is_srr] <- unname(run_to_sample[tip_new[is_srr]])

  tree$tip.label <- make.unique(tip_new, sep = "_dup")
  tree
}

map_mitotypes_to_sample_only <- function(mitotypes_tbl, run_to_sample) {
  mitotypes_tbl %>%
    mutate(sample = if_else(sample %in% names(run_to_sample),
                            unname(run_to_sample[sample]),
                            sample)) %>%
    mutate(sample = make.unique(sample, sep = "_dup")) %>%
    distinct(sample, .keep_all = TRUE)
}

# prune until prop_diff has no NA off-diagonal (needed for NJ)
prune_to_complete_dist <- function(prop_diff, compared, min_keep = 3) {
  keep <- colnames(prop_diff)
  pd <- prop_diff
  cp <- compared

  repeat {
    offdiag_na <- is.na(pd)
    diag(offdiag_na) <- FALSE
    if (!any(offdiag_na)) break

    na_per_sample <- colSums(offdiag_na) + rowSums(offdiag_na)
    worst <- names(which.max(na_per_sample))

    keep <- setdiff(keep, worst)
    if (length(keep) < min_keep) {
      break
    }
    pd <- pd[keep, keep, drop = FALSE]
    cp <- cp[keep, keep, drop = FALSE]
  }

  list(keep = keep, prop_diff = pd, compared = cp)
}

# ============================================================
# 2) READ METADATA ONCE
# ============================================================

samples_tbl <- read_csv(samples_csv, show_col_types = FALSE)

# SRA table for D. obtusa
sra_map_all <- read_tsv(sra_tsv, show_col_types = FALSE) %>%
  select(Run, ScientificName, SampleName) %>%
  filter(!is.na(Run), !is.na(SampleName), SampleName != "") %>%
  filter(ScientificName == "Daphnia obtusa")

run_to_sample <- setNames(sra_map_all$SampleName, sra_map_all$Run)

# Parse SRA location from SampleName
sra_loc_tbl <- sra_map_all %>%
  mutate(
    sra_location = parse_sra_location_from_samplename(SampleName),
    sra_location_canon = canon_loc(sra_location)
  ) %>%
  filter(!is.na(sra_location), sra_location != "")

# Local metadata with clone column, define location as "prefix before _"
if (!is.na(metadata_with_clone_csv) && file.exists(metadata_with_clone_csv)) {
  tmp <- read_csv(metadata_with_clone_csv, show_col_types = FALSE)
  if (!all(c("Well", "clone") %in% names(tmp))) {
    stop("metadata_with_clone_csv must contain columns: Well, clone")
  }
  local_loc_tbl <- tmp %>%
    transmute(
      Well = trimws(as.character(Well)),
      clone = trimws(as.character(clone)),
      loc_prefix = clone_prefix_before_underscore(clone),
      loc_prefix_canon = canon_loc(loc_prefix)
    ) %>%
    filter(!is.na(loc_prefix), loc_prefix != "", !is.na(Well), Well != "")
} else {
  stop("metadata_with_clone_csv not found: ", metadata_with_clone_csv)
}

# ============================================================
# 3) MULTI-LOCATION (combined) CASE-INSENSITIVE: COLLECT TARGET SAMPLES
# ============================================================

locs_c <- sort(unique(canon_loc(LOCATIONS)))
locs_c <- sort(unique(canon_loc(LOCATIONS)))

# short + stable-ish label (doesn't exceed filename limits)
LOCATION_LABEL <- paste0(
  "locset_n", length(locs_c),
  "__", paste(head(locs_c, 3), collapse = "-"),
  if (length(locs_c) > 3) "__etc" else ""
)

out_dir <- file.path(BASE_OUT, LOCATION_LABEL)

# IMPORTANT: don't hide warnings, and error if it fails
ok <- dir.create(out_dir, recursive = TRUE, showWarnings = TRUE)
if (!dir.exists(out_dir)) {
  stop("Failed to create out_dir: ", out_dir,
       "\n(Usually this is a path length or permissions issue.)")
}

message("Output dir = ", out_dir)
message("\n==============================")
message("LOCATIONS (canonical) = ", paste(locs_c, collapse = ", "))
message("Output dir           = ", out_dir)
message("==============================")

local_loc_ids <- local_loc_tbl %>%
  filter(loc_prefix_canon %in% locs_c) %>%
  pull(Well) %>%
  unique()

sra_runs_location <- sra_loc_tbl %>%
  filter(sra_location_canon %in% locs_c) %>%
  pull(Run) %>%
  unique()

matched_local_locs <- sort(unique(local_loc_tbl$loc_prefix[local_loc_tbl$loc_prefix_canon %in% locs_c]))
matched_sra_locs   <- sort(unique(sra_loc_tbl$sra_location[sra_loc_tbl$sra_location_canon %in% locs_c]))
message("Matched local loc_prefix values: ", ifelse(length(matched_local_locs) == 0, "(none)", paste(matched_local_locs, collapse = ", ")))
message("Matched SRA locations:          ", ifelse(length(matched_sra_locs) == 0, "(none)", paste(matched_sra_locs, collapse = ", ")))

target_samples <- unique(c(local_loc_ids, sra_runs_location))

message("Target samples (selected locations): ", length(target_samples),
        " [local wells=", length(local_loc_ids), ", SRA runs=", length(sra_runs_location), "]")

if (length(target_samples) < 3) {
  stop("Too few target samples before GDS intersection. Check LOCATIONS vs metadata.")
}

# ============================================================
# 4) OPEN GDS, FILTER, LOAD GENOTYPES
# ============================================================

gds <- seqOpen(gds_path)
on.exit(seqClose(gds), add = TRUE)

seqResetFilter(gds)
gds_samples <- seqGetData(gds, "sample.id")
target_samples_in_gds <- intersect(target_samples, gds_samples)

if (length(target_samples_in_gds) < 3) {
  stop("Too few target samples in GDS for selected LOCATIONS (", length(target_samples_in_gds), ")")
}

seqSetFilter(gds, sample.id = target_samples_in_gds, verbose = FALSE)

# Per-sample missingness
miss_sample <- seqMissing(gds, per.variant = FALSE)
sample_ids  <- seqGetData(gds, "sample.id")
keep_samp1  <- sample_ids[miss_sample < MISS_SAMPLE_MAX]
if (length(keep_samp1) < 3) stop("Too few samples after missingness filter")
seqSetFilter(gds, sample.id = keep_samp1, verbose = FALSE)

# Per-variant missingness
miss_var <- seqMissing(gds, per.variant = TRUE)
var_ids  <- seqGetData(gds, "variant.id")
keep_var <- var_ids[miss_var < MISS_SITE_MAX]
if (length(keep_var) < 50) stop("Too few variants after variant filter")
seqSetFilter(gds, sample.id = keep_samp1, variant.id = keep_var, verbose = FALSE)

# Mean DP per sample
DP_raw <- seqGetData(gds, "annotation/format/DP")
sample_ids2 <- seqGetData(gds, "sample.id")
mean_dp <- mean_dp_per_sample(DP_raw, sample_ids2)

keep_samp2 <- sample_ids2[!is.na(mean_dp) & mean_dp >= MIN_MEAN_DP]
if (length(keep_samp2) < 3) stop("Too few samples after DP filter")
seqSetFilter(gds, sample.id = keep_samp2, variant.id = keep_var, verbose = FALSE)

message("Kept samples:  ", length(seqGetData(gds, "sample.id")))
message("Kept variants: ", length(seqGetData(gds, "variant.id")))

# Load genotype (variants x samples)
GT_raw <- seqGetData(gds, "genotype")
stopifnot(length(dim(GT_raw)) == 3L)

d <- dim(GT_raw)
nsamp <- length(seqGetData(gds, "sample.id"))

if (d[2] == nsamp) {
  GT <- t(GT_raw[1, , ])
} else if (d[3] == nsamp) {
  GT <- GT_raw[1, , ]
} else {
  stop("Cannot infer GT dim order; dim=", paste(d, collapse = "x"))
}

GT <- as.matrix(GT)
colnames(GT) <- seqGetData(gds, "sample.id")

# Collapse het/hom alt into 1; keep missing as-is (handled in pairwise_diffs)
GT[!is.na(GT) & GT >= 0L & GT > 1] <- 1L

# ============================================================
# 5) PAIRWISE DIFFS + OVERLAP FILTER (prevents empty tree)
# ============================================================

pw0 <- pairwise_diffs(GT)
compared0 <- pw0$compared

ok_pairs_per_sample <- rowSums(compared0 >= MIN_PAIRWISE_COMPARED, na.rm = TRUE)
keep_overlap <- names(ok_pairs_per_sample)[ok_pairs_per_sample >= MIN_OK_PAIRS_PER_SAMPLE]

message("Samples with >= ", MIN_OK_PAIRS_PER_SAMPLE,
        " partners at overlap >= ", MIN_PAIRWISE_COMPARED, ": ", length(keep_overlap))

if (length(keep_overlap) < MIN_SAMPLES_FOR_TREE) {
  stop("After overlap filtering, too few samples remain (", length(keep_overlap),
       "). Try lowering MIN_PAIRWISE_COMPARED or MIN_OK_PAIRS_PER_SAMPLE.")
}

GT <- GT[, keep_overlap, drop = FALSE]

pw <- pairwise_diffs(GT)
diff_counts <- pw$diff_counts
compared    <- pw$compared

# prop_diff with NA where compared==0
prop_diff <- matrix(NA_real_, nrow(diff_counts), ncol(diff_counts),
                    dimnames = dimnames(diff_counts))
idx_ok <- compared > 0
prop_diff[idx_ok] <- diff_counts[idx_ok] / compared[idx_ok]

# Ensure we have a complete distance matrix for NJ (no NA off-diagonal)
pr <- prune_to_complete_dist(prop_diff, compared, min_keep = MIN_SAMPLES_FOR_TREE)
keep_complete <- pr$keep
prop_diff_complete <- pr$prop_diff
compared_complete  <- pr$compared

if (length(keep_complete) < MIN_SAMPLES_FOR_TREE) {
  stop("Cannot form a complete distance matrix for NJ after pruning. ",
       "This usually means too little overlap among samples for these locations. ",
       "Try lowering MIN_PAIRWISE_COMPARED.")
}

# Subset GT + matrices to the complete set used for NJ + mitotypes
GT <- GT[, keep_complete, drop = FALSE]
diff_counts <- diff_counts[keep_complete, keep_complete, drop = FALSE]
compared    <- compared[keep_complete, keep_complete, drop = FALSE]
prop_diff   <- prop_diff_complete

message("Final samples used for distances/tree: ", length(keep_complete))

# Save matrices
saveRDS(list(diff_counts = diff_counts, compared = compared, prop_diff = prop_diff),
        file = file.path(out_dir, "allsite_pairwise_matrices.rds"))

# Long tables
DT_prop <- as.data.table(as.table(prop_diff));    setnames(DT_prop, c("sampleA","sampleB","prop_diff"))
DT_diff <- as.data.table(as.table(diff_counts));  setnames(DT_diff, c("sampleA","sampleB","diffs"))
DT_comp <- as.data.table(as.table(compared));     setnames(DT_comp, c("sampleA","sampleB","compared_sites"))

pairwise_long <- merge(DT_prop, DT_diff, by = c("sampleA","sampleB"))
pairwise_long <- merge(pairwise_long, DT_comp, by = c("sampleA","sampleB"))
pairwise_long <- pairwise_long[sampleA != sampleB]
pairwise_long[, Similarity := 1 - prop_diff]
fwrite(pairwise_long, file.path(out_dir, "allsite_pairwise_long.tsv"), sep = "\t")

# Diff-rate plots (skip if no finite diffs)
finite_diffs <- pairwise_long[is.finite(diffs)]
if (nrow(finite_diffs) > 0) {
  p1 <- ggplot(finite_diffs, aes(x = diffs)) +
    geom_histogram(bins = 100) + theme_bw() + xlim(0, 100) +
    labs(x = "diffs", y = "Frequency", title = paste0("Pairwise diffs (<=100), ", LOCATION_LABEL))
  ggsave(file.path(out_dir, "diff_rates.png"), p1, width = 7, height = 6, dpi = 300)

  p2 <- ggplot(finite_diffs, aes(x = diffs)) +
    geom_histogram(bins = 100) + theme_bw() +
    labs(x = "diffs", y = "Frequency", title = paste0("Pairwise diffs (all), ", LOCATION_LABEL))
  ggsave(file.path(out_dir, "diff_rates_all.png"), p2, width = 7, height = 6, dpi = 300)

  p3 <- ggplot(finite_diffs, aes(x = diffs)) +
    geom_histogram(bins = 300) + theme_bw() + xlim(0, 300) +
    labs(x = "diffs", y = "Frequency", title = paste0("Pairwise diffs (<=300), ", LOCATION_LABEL))
  ggsave(file.path(out_dir, "diff_rates_300.png"), p3, width = 7, height = 6, dpi = 300)
} else {
  message("WARN: no finite diffs to plot.")
}

# ============================================================
# 6) ASSIGN MITOTYPES (< DIFF_THRESHOLD) + SAVE
# ============================================================

group_letters <- assign_groups_from_diffs(diff_counts, compared, threshold = DIFF_THRESHOLD)

mitotypes_tbl <- tibble(sample = names(group_letters), mitotype = unname(group_letters)) %>%
  arrange(mitotype, sample)

# If somehow empty, stop loudly (better than silently writing empty outputs)
if (nrow(mitotypes_tbl) == 0) {
  stop("Mitotypes table is empty. This indicates the graph had no edges (no pairs with diffs < threshold). ",
       "Try increasing DIFF_THRESHOLD or verify prop_diff is reasonable.")
}

write_csv(mitotypes_tbl, file.path(out_dir, "eudobtusa_mito_types.csv"))

pairwise_long2 <- as_tibble(pairwise_long) %>%
  mutate(
    mitoA = mitotypes_tbl$mitotype[match(sampleA, mitotypes_tbl$sample)],
    mitoB = mitotypes_tbl$mitotype[match(sampleB, mitotypes_tbl$sample)]
  )
write_csv(pairwise_long2, file.path(out_dir, "allsite_pairwise_long_with_mitotypes.csv"))

# ============================================================
# 7) MITOTYPE PAIRWISE COMPARISONS
# ============================================================

pairwise_tri <- pairwise_long2 %>%
  filter(sampleA < sampleB) %>%
  filter(!is.na(mitoA), !is.na(mitoB)) %>%
  mutate(comp_type = if_else(mitoA == mitoB, "Within mitotype", "Between mitotypes"))

if (nrow(pairwise_tri) > 0) {
  p_within_between <- ggplot(pairwise_tri, aes(x = comp_type, y = prop_diff)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15, height = 0, alpha = 0.4) +
    theme_bw() +
    labs(title = paste0("Within vs Between mitotype distances, ", LOCATION_LABEL),
         x = "", y = "Proportion different")
  ggsave(file.path(out_dir, "mitotype_within_vs_between_box.png"),
         p_within_between, width = 7, height = 5, dpi = 300)

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
    labs(title = paste0("Pairwise distances by mitotype pair, ", LOCATION_LABEL),
         x = "Mitotype pair", y = "Proportion different")
  ggsave(file.path(out_dir, "mitotype_pairwise_box.png"),
         p_pairs, width = 10, height = 6, dpi = 300)

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
} else {
  message("WARN: No comparable mitotype pairs after filtering.")
}

# ============================================================
# 8) NJ TREE + ROOT + REPLACE SRR WITH SampleName ONLY
# ============================================================

# prop_diff should have NO NA off-diagonal here
if (any(is.na(prop_diff[upper.tri(prop_diff)]))) {
  stop("Distance matrix still has NA values; cannot build NJ tree.")
}

dist_matrix <- as.dist(prop_diff)
tree <- nj(dist_matrix)

if (ROOT_OUTGROUP %in% tree$tip.label) {
  tree_rooted <- root(tree, outgroup = ROOT_OUTGROUP, resolve.root = TRUE)
} else {
  tree_rooted <- tree
}
tree_rooted$edge.length[tree_rooted$edge.length < 0] <- 0

tree_rooted2 <- make_tree_tiplabels_sample_only(tree_rooted, run_to_sample)
write.tree(tree_rooted2, file = file.path(out_dir, "eudobtusa_tree_nwk.nwk"))

# color by mitotype in SampleName space
mito_sample_only <- map_mitotypes_to_sample_only(mitotypes_tbl, run_to_sample)
mito_map <- setNames(mito_sample_only$mitotype, mito_sample_only$sample)

tip_mito <- unname(mito_map[tree_rooted2$tip.label])
mito_levels <- sort(unique(na.omit(tip_mito)))
mito_cols <- setNames(rainbow(length(mito_levels)), mito_levels)

tip_colors <- mito_cols[tip_mito]
tip_colors[is.na(tip_colors)] <- "grey50"

xmax <- as.numeric(quantile(node.depth.edgelength(tree_rooted2), 0.98))







clone_tbl <- local_loc_tbl %>%
  distinct(Well, clone) %>%
  mutate(
    Well  = trimws(as.character(Well)),
    clone = trimws(as.character(clone))
  ) %>%
  filter(!is.na(Well), Well != "", !is.na(clone), clone != "") %>%
  mutate(
    well_token   = sub("^.*_", "", Well),
    prefix_token = if_else(grepl("_", Well), sub("_[^_]+$", "", Well), NA_character_)
  )

tip_ids <- tree_rooted2$tip.label

tip_tbl <- tibble(tip = tip_ids) %>%
  mutate(
    well_token   = sub("^.*_", "", tip),
    prefix_token = if_else(grepl("_", tip), sub("_[^_]+$", "", tip), NA_character_),
    want_clone   = grepl("^(RobertUK|Gilmer|Rockpool)", tip)
  )

tip_labeled <- tip_tbl %>%
  left_join(
    clone_tbl %>% select(prefix_token, well_token, clone),
    by = c("prefix_token", "well_token")
  ) %>%
  mutate(clone2 = clone) %>%
  select(-clone) %>%
  left_join(
    clone_tbl %>% transmute(tip = Well, clone_exact = clone),
    by = "tip"
  ) %>%
  mutate(clone_final = coalesce(clone2, clone_exact))

display_label <- tip_labeled$tip
use_clone <- tip_labeled$want_clone & !is.na(tip_labeled$clone_final)
display_label[use_clone] <- tip_labeled$clone_final[use_clone]
display_label <- make.unique(display_label, sep = "_dup")

# Debug: save the mapping you’re expecting to see on the plot
write_csv(
  tip_labeled %>% transmute(tip, display_label, clone_final),
  file.path(out_dir, "tip_label_mapping.csv")
)

# ---- IMPORTANT: give the plot extra radial space so labels are not clipped ----
depths <- node.depth.edgelength(tree_rooted2)
xmax_tip <- max(depths, na.rm = TRUE)
xmax_plot <- xmax_tip * 1.25   # extra room for text
label_offset <- xmax_tip * 0.03  # push labels out a bit

# Write NEW filenames so you can tell you’re looking at the updated output
pdf(file.path(out_dir, "eudobtusa_groups_circ_CLONELABELS.pdf"), width = 30, height = 16)
plot.phylo(
  tree_rooted2,
  type = "fan",
  cex = 0.7,
  no.margin = TRUE,
  tip.color = tip_colors,
  x.lim = c(0, xmax_plot),
  show.tip.label = FALSE,
  main = paste0("Rooted NJ Tree (mito, all sites), ", LOCATION_LABEL)
)

tiplabels(
  text   = display_label,
  tip    = seq_along(display_label),
  frame  = "none",
  cex    = 0.7,
  col    = tip_colors,
  adj    = c(0, 0.5),
  offset = label_offset
)

legend("topleft", legend = mito_levels, col = mito_cols[mito_levels], pch = 19,
       pt.cex = 1.3, cex = 0.9, bty = "n", title = "Mitotype")
dev.off()

png(file.path(out_dir, "eudobtusa_groups_circ_CLONELABELS.png"), width = 3500, height = 2750, res = 300)
plot.phylo(
  tree_rooted2,
  type = "fan",
  cex = 0.7,
  no.margin = TRUE,
  tip.color = tip_colors,
  show.tip.label = FALSE,
  main = paste0("Rooted NJ Tree (mito, all sites), ", LOCATION_LABEL)
)

tiplabels(
  text   = display_label,
  tip    = seq_along(display_label),
  frame  = "none",
  cex    = 0.7,
  col    = tip_colors,
  adj    = c(0, 0.5),
  offset = label_offset
)

legend("topleft", legend = mito_levels, col = mito_cols[mito_levels], pch = 19,
       pt.cex = 1.3, cex = 0.9, bty = "n", title = "Mitotype")
dev.off()


# ============================================================
# 9) SAVE KEPT LISTS
# ============================================================

writeLines(seqGetData(gds, "sample.id"), file.path(out_dir, "kept_samples.txt"))
writeLines(as.character(seqGetData(gds, "variant.id")), file.path(out_dir, "kept_variants.txt"))

# ============================================================
# 10) REPRESENTATIVE PER MITOTYPE (highest mean DP)
# ============================================================

DP_raw2 <- seqGetData(gds, "annotation/format/DP")
kept_ids_now <- seqGetData(gds, "sample.id")
mean_dp_now <- mean_dp_per_sample(DP_raw2, kept_ids_now)

dp_tbl <- tibble(sample = kept_ids_now, meanDP = as.numeric(mean_dp_now))

dp_mito <- dp_tbl %>%
  left_join(mitotypes_tbl, by = "sample") %>%
  filter(!is.na(mitotype))

if (nrow(dp_mito) == 0) {
  message("WARN: no mitotype-assigned samples for representatives")
} else {
  rep_per_mitotype <- dp_mito %>%
    arrange(mitotype, desc(meanDP), sample) %>%
    group_by(mitotype) %>%
    slice(1) %>%
    ungroup()

  write_csv(rep_per_mitotype, file.path(out_dir, "representative_per_mitotype_highestDP.csv"))
  writeLines(rep_per_mitotype$sample, file.path(out_dir, "representative_per_mitotype_highestDP.txt"))

  rep_per_mitotype_display <- rep_per_mitotype %>%
    mutate(display_name = if_else(sample %in% names(run_to_sample),
                                  unname(run_to_sample[sample]),
                                  sample),
           Location = LOCATION_LABEL)

  write_csv(rep_per_mitotype_display,
            file.path(out_dir, "representative_per_mitotype_highestDP_withSampleName.csv"))
}

message("DONE locations: ", LOCATION_LABEL)












# ============================================================
# 10) REPRESENTATIVES:
#     - top 3 per mitotype (or fewer if <3 exist)
#     - AND at least 1 per LOCAL clone-location (prefix before "_" in clone)
#     - AND at least 1 per SRA location (middle token in SampleName)
# ============================================================

TOP_N_PER_MITOTYPE <- 3

DP_raw2 <- seqGetData(gds, "annotation/format/DP")
kept_ids_now <- seqGetData(gds, "sample.id")
mean_dp_now <- mean_dp_per_sample(DP_raw2, kept_ids_now)

dp_tbl <- tibble(sample = kept_ids_now, meanDP = as.numeric(mean_dp_now))

dp_mito <- dp_tbl %>%
  left_join(mitotypes_tbl, by = "sample") %>%
  filter(!is.na(mitotype))

if (nrow(dp_mito) == 0) {
  message("WARN: no mitotype-assigned samples for representatives")
} else {

  # ---- Local Well -> clone map ----
  clone_tbl <- local_loc_tbl %>%
    distinct(Well, clone) %>%
    mutate(
      Well  = trimws(as.character(Well)),
      clone = trimws(as.character(clone))
    ) %>%
    filter(!is.na(Well), Well != "", !is.na(clone), clone != "")

  clone_map <- setNames(clone_tbl$clone, clone_tbl$Well)

  # ---- annotate dp_mito with:
  #      - display_name (SRA SampleName if SRR, else sample)
  #      - sra_location (parsed from SampleName)
  #      - clone + clone_loc (prefix before "_" in clone)
  dp_mito2 <- dp_mito %>%
    mutate(
      display_name = if_else(sample %in% names(run_to_sample),
                             unname(run_to_sample[sample]),
                             sample),

      # SRA location from SampleName (only makes sense for the SRA display_name format)
      sra_location = parse_sra_location_from_samplename(display_name),
      sra_location_canon = canon_loc(sra_location),

      # Local clone lookup (only for local samples that exist in clone_map keys)
      clone = unname(clone_map[sample]),
      clone_loc = if_else(!is.na(clone) & grepl("_", clone),
                          sub("_.*$", "", clone),
                          NA_character_),
      clone_loc_canon = canon_loc(clone_loc)
    )

  # ------------------------------------------------------------
  # A) top N per mitotype
  # ------------------------------------------------------------
  reps_topN <- dp_mito2 %>%
    arrange(mitotype, desc(meanDP), sample) %>%
    group_by(mitotype) %>%
    slice_head(n = TOP_N_PER_MITOTYPE) %>%
    ungroup()

  # ------------------------------------------------------------
  # B) at least 1 per LOCAL clone-location
  # ------------------------------------------------------------
  reps_one_per_local_loc <- dp_mito2 %>%
    filter(!is.na(clone_loc_canon), clone_loc_canon != "") %>%
    arrange(clone_loc_canon, desc(meanDP), sample) %>%
    group_by(clone_loc_canon) %>%
    slice_head(n = 1) %>%
    ungroup()

  # ------------------------------------------------------------
  # C) at least 1 per SRA location
  # ------------------------------------------------------------
  reps_one_per_sra_loc <- dp_mito2 %>%
    filter(!is.na(sra_location_canon), sra_location_canon != "") %>%
    arrange(sra_location_canon, desc(meanDP), sample) %>%
    group_by(sra_location_canon) %>%
    slice_head(n = 1) %>%
    ungroup()

  # ------------------------------------------------------------
  # D) union + deduplicate
  # ------------------------------------------------------------
  reps_final <- bind_rows(reps_topN, reps_one_per_local_loc, reps_one_per_sra_loc) %>%
    distinct(sample, .keep_all = TRUE) %>%
    arrange(mitotype, desc(meanDP), sample)

  # ---- Save outputs ----
  write_csv(reps_final, file.path(out_dir, "representatives_top3_per_mitotype_plus1_per_localLoc_plus1_per_sraLoc.csv"))
  writeLines(reps_final$sample, file.path(out_dir, "representatives_top3_per_mitotype_plus1_per_localLoc_plus1_per_sraLoc.txt"))

  reps_final_display <- reps_final %>%
    mutate(LocationSet = LOCATION_LABEL)

  write_csv(
    reps_final_display,
    file.path(out_dir, "representatives_top3_per_mitotype_plus1_per_localLoc_plus1_per_sraLoc_withDisplayName.csv")
  )

  # Optional summaries to sanity-check coverage
  write_csv(reps_final %>% filter(!is.na(clone_loc_canon), clone_loc_canon != "") %>% count(clone_loc_canon),
            file.path(out_dir, "representatives_covered_local_locations.csv"))
  write_csv(reps_final %>% filter(!is.na(sra_location_canon), sra_location_canon != "") %>% count(sra_location_canon),
            file.path(out_dir, "representatives_covered_sra_locations.csv"))
  write_csv(reps_final %>% count(mitotype),
            file.path(out_dir, "representatives_summary_by_mitotype.csv"))
}

message("DONE locations: ", LOCATION_LABEL)








