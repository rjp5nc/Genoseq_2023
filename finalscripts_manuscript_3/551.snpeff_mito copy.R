#ijob -A berglandlab -c10 -p standard --mem=40G
#module load gcc/11.4.0  openmpi/4.1.4 icu R/4.3.1;R
#R

suppressPackageStartupMessages({
  library(vcfR)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(readr)
})

# ============================================================
# USER SETTINGS (EDIT THESE)
# ============================================================
VCF_IN   <- "/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/snpeff/output.ann.vcf"
MITO_TSV <- "/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/mito_types.csv"
OUT_DIR  <- "/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/snpeff"
PREFIX   <- "biallelic"   # filename prefix for all outputs

# Optional: if your mitotype file does not already have headers "sample,mitotype",
# set this TRUE to force the first two columns to be sample/mitotype
FORCE_TWO_COLS_SAMPLE_MITOTYPE <- TRUE

# Optional: drop tiny totals to reduce noise in proportion plots
MIN_TOTAL_WITHIN <- 1
MIN_TOTAL_PAIR   <- 1

# ============================================================
# OUTPUT FILES
# ============================================================
stopifnot(dir.exists(OUT_DIR))
OUT_WITHIN_TSV <- file.path(OUT_DIR, paste0(PREFIX, ".within.props.tsv"))
OUT_WITHIN_PNG <- file.path(OUT_DIR, paste0(PREFIX, ".within.props.png"))

OUT_PAIR_TSV   <- file.path(OUT_DIR, paste0(PREFIX, ".pairwise.between.props.tsv"))
OUT_PAIR_PNG   <- file.path(OUT_DIR, paste0(PREFIX, ".pairwise.between.props.png"))

OUT_BAR_TSV    <- file.path(OUT_DIR, paste0(PREFIX, ".within.counts.tsv"))
OUT_BAR_PNG    <- file.path(OUT_DIR, paste0(PREFIX, ".within.counts.png"))

OUT_LOG        <- file.path(OUT_DIR, paste0(PREFIX, ".log.txt"))

log_msg <- function(...) {
  msg <- paste0(...)
  message(msg)
  cat(msg, "\n", file = OUT_LOG, append = TRUE)
}

# ============================================================
# LOAD MITOTYPE MAP
# ============================================================
mito_df_raw <- suppressMessages(readr::read_delim(
  MITO_TSV,
  delim = ifelse(grepl("\\.csv$", MITO_TSV, ignore.case = TRUE), ",", "\t"),
  show_col_types = FALSE
))

if (FORCE_TWO_COLS_SAMPLE_MITOTYPE) {
  if (ncol(mito_df_raw) < 2) stop("MITO_TSV has <2 columns.")
  mito_df_raw <- mito_df_raw[, 1:2]
  colnames(mito_df_raw) <- c("sample", "mitotype")
}

# If already named correctly, keep; else try a few common patterns
if (all(c("sample", "mitotype") %in% names(mito_df_raw))) {
  mito_df <- mito_df_raw %>% select(sample, mitotype)
} else if (all(c("Well", "clone") %in% names(mito_df_raw))) {
  mito_df <- mito_df_raw %>% transmute(sample = clone, mitotype = Well)
} else if (all(c("sample", "haplotype") %in% names(mito_df_raw))) {
  mito_df <- mito_df_raw %>% transmute(sample = sample, mitotype = haplotype)
} else {
  stop("Couldn't infer mitotype mapping columns. Make a 2-col file: sample,mitotype")
}

mito_df <- mito_df %>%
  filter(!is.na(sample), !is.na(mitotype)) %>%
  mutate(sample = as.character(sample),
         mitotype = as.character(mitotype)) %>%
  distinct()

log_msg("Loaded mitotype map: ", nrow(mito_df), " rows; ", length(unique(mito_df$mitotype)), " mitotypes.")

# ============================================================
# READ VCF
# ============================================================
vcf <- vcfR::read.vcfR(VCF_IN, verbose = FALSE)
gt  <- vcfR::extract.gt(vcf, element = "GT") # variants x samples
fix <- as.data.frame(vcf@fix, stringsAsFactors = FALSE)

# Keep samples with mitotype labels
keep_samples <- intersect(colnames(gt), mito_df$sample)
if (length(keep_samples) < 2) stop("Too few VCF samples match mitotype map. Check sample naming.")
gt <- gt[, keep_samples, drop = FALSE]
mito_df <- mito_df %>% filter(sample %in% keep_samples)

log_msg("VCF variants: ", nrow(fix), " ; samples matched to mitotypes: ", length(keep_samples))

# ============================================================
# KEEP BIALLELIC SITES (ALT has no comma) AND SPLIT MULTI-ANN?
# Here we treat "biallelic record" as ALT field with exactly one allele
# ============================================================
bial_mask <- !is.na(fix$ALT) & !str_detect(fix$ALT, ",")
fix_bi <- fix[bial_mask, , drop = FALSE]
gt_bi  <- gt[bial_mask, , drop = FALSE]
log_msg("Biallelic records kept: ", nrow(fix_bi))

# ============================================================
# PARSE snpEff ANN -> gene + effect (syn/nonsyn)
# NOTE: Return MUST always be character(2) to satisfy vapply
# We prefer nonsyn if present anywhere in ANN entries
# ============================================================
parse_ann <- function(info) {
  if (is.na(info) || !str_detect(info, "ANN=")) return(c(NA_character_, NA_character_))

  ann <- str_replace(info, "^.*ANN=", "")
  ann <- str_split(ann, ";", n = 2, simplify = TRUE)[1]  # up to first ';'
  entries <- str_split(ann, ",")[[1]]

  best_gene <- NA_character_
  best_eff  <- NA_character_

  for (e in entries) {
    f <- str_split(e, "\\|", simplify = TRUE)
    if (ncol(f) < 5) next
    eff  <- f[1, 2]
    gene <- f[1, 4]
    if (is.na(gene) || gene == "") next

    # treat combined annotations like "missense_variant&splice_region_variant"
    if (str_detect(eff, "missense_variant")) {
      return(c(gene, "nonsyn"))
    }
    if (str_detect(eff, "synonymous_variant") && is.na(best_eff)) {
      best_gene <- gene
      best_eff  <- "syn"
      # keep looking in case a missense shows up later
    }
  }

  c(best_gene, best_eff)
}

anno_mat <- t(vapply(fix_bi$INFO, parse_ann, character(2)))
fix_bi$gene   <- anno_mat[, 1]
fix_bi$effect <- anno_mat[, 2]

# Keep only syn/nonsyn with gene name
keep_eff <- !is.na(fix_bi$gene) & fix_bi$effect %in% c("syn", "nonsyn")
fix2 <- fix_bi[keep_eff, , drop = FALSE]
gt2  <- gt_bi[keep_eff, , drop = FALSE]
log_msg("Annotated syn/nonsyn biallelic records kept: ", nrow(fix2))
log_msg("Effect table:\n", paste(capture.output(print(table(fix2$effect))), collapse = "\n"))
log_msg("Unique genes: ", length(unique(fix2$gene)))

# ============================================================
# GENOTYPE -> ALLELE (haploid-ish: take first allele before / or |)
# ============================================================
gt_to_allele <- function(x) {
  if (is.na(x) || x %in% c(".", "./.", ".|.")) return(NA_integer_)
  a <- str_split(x, "[/|]", simplify = TRUE)[1]
  suppressWarnings(as.integer(a))
}


first_allele <- function(x) {
  x <- sub("[/|].*$", "", x)          # keep first allele
  x[x %in% c(".", "./.", ".|.")] <- NA
  suppressWarnings(as.integer(x))
}

allele_mat <- matrix(first_allele(gt2),
                     nrow = nrow(gt2),
                     ncol = ncol(gt2),
                     dimnames = dimnames(gt2))


log_msg("Non-missing allele calls: ", sum(!is.na(allele_mat)))

# ============================================================
# GROUPS BY MITOTYPE
# ============================================================
mito_levels     <- sort(unique(mito_df$mitotype))
mito_to_samples <- split(mito_df$sample, mito_df$mitotype)

# Fixed allele per mitotype and polymorphic flag
fixed_by_mito <- matrix(NA_integer_, nrow=nrow(allele_mat), ncol=length(mito_levels),
                        dimnames=list(NULL, mito_levels))
poly_by_mito  <- matrix(FALSE, nrow=nrow(allele_mat), ncol=length(mito_levels),
                        dimnames=list(NULL, mito_levels))

for (m in mito_levels) {
  s <- mito_to_samples[[m]]
  sub <- allele_mat[, s, drop = FALSE]
  aset <- apply(sub, 1, function(v) sort(unique(v[!is.na(v)])))
  poly_by_mito[, m]  <- vapply(aset, function(a) length(a) > 1, logical(1))
  fixed_by_mito[, m] <- vapply(aset, function(a) if (length(a)==1) a[1] else NA_integer_, integer(1))
}

# Variant annotation aligned with allele_mat rows
var_anno <- tibble(
  variant_index = seq_len(nrow(fix2)),
  gene = fix2$gene,
  effect = fix2$effect
)
stopifnot(nrow(var_anno) == nrow(allele_mat))

# ============================================================
# (A) WITHIN-MITOTYPE: proportions syn vs nonsyn by gene
# ============================================================
within_long <- as.data.frame(poly_by_mito) %>%
  mutate(variant_index = seq_len(nrow(poly_by_mito))) %>%
  pivot_longer(-variant_index, names_to="mitotype", values_to="is_poly") %>%
  filter(is_poly) %>%
  select(variant_index, mitotype) %>%
  left_join(var_anno, by="variant_index") %>%
  filter(!is.na(gene), effect %in% c("syn","nonsyn"))

within_props <- within_long %>%
  count(mitotype, gene, effect, name="n") %>%
  pivot_wider(names_from = effect, values_from = n, values_fill = 0) %>%
  mutate(
    total = syn + nonsyn,
    prop_nonsyn = ifelse(total > 0, nonsyn / total, NA_real_),
    prop_syn    = ifelse(total > 0, syn    / total, NA_real_)
  ) %>%
  filter(total >= MIN_TOTAL_WITHIN) %>%
  arrange(gene, mitotype)

write_tsv(within_props, OUT_WITHIN_TSV)
log_msg("Wrote: ", OUT_WITHIN_TSV)

p_within <- ggplot(within_props, aes(x = prop_nonsyn, y = prop_syn, color = mitotype)) +
  geom_point(size = 2, alpha = 0.9) +
  coord_equal(xlim = c(0,1), ylim = c(0,1)) +
  facet_wrap(~ gene) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  labs(
    x = "Proportion nonsynonymous (within mitotype)",
    y = "Proportion synonymous (within mitotype)",
    color = "Mitotype",
    title = "Within-mitotype syn vs nonsyn proportions by gene (biallelic)"
  )

ggsave(OUT_WITHIN_PNG, p_within, width = 14, height = 10, dpi = 300)
log_msg("Wrote: ", OUT_WITHIN_PNG)


OUT_LOC_WITHIN_TSV <- file.path(OUT_DIR, paste0(PREFIX, ".within.location_mitotype.PNPS_by_gene.tsv"))
OUT_LOC_WITHIN_PNG <- file.path(OUT_DIR, paste0(PREFIX, ".within.location_mitotype.PNPS_by_gene.png"))

# ---- sample -> accuratelocation (standardize metadata clone key)
meta_loc <- metadata %>%
  transmute(
    sample = str_remove(as.character(clone), "_clone$"),
    accuratelocation = na_if(as.character(accuratelocation), "")
  ) %>%
  distinct()

# join location onto mito_df (which is the VCF sample list + mitotype)
sample_map <- mito_df %>%
  mutate(sample = as.character(sample),
         mitotype = as.character(mitotype)) %>%
  left_join(meta_loc, by = "sample") %>%
  filter(!is.na(accuratelocation), accuratelocation != "")

# subset allele_mat columns to samples that have locations
keep_loc_samples <- intersect(colnames(allele_mat), sample_map$sample)
allele_loc <- allele_mat[, keep_loc_samples, drop = FALSE]
sample_map <- sample_map %>% filter(sample %in% keep_loc_samples)

log_msg("Samples with accuratelocation: ", nrow(sample_map),
        " across ", n_distinct(sample_map$accuratelocation), " locations.")

# ---- helper: polymorphic within a set of samples
is_poly_within <- function(allele_mat, samples) {
  if (length(samples) < 2) return(rep(FALSE, nrow(allele_mat)))
  sub <- allele_mat[, samples, drop = FALSE]
  apply(sub, 1, function(v) {
    vv <- v[!is.na(v)]
    length(unique(vv)) > 1
  })
}

# ---- compute poly flags per (location, mitotype)
loc_levels  <- sort(unique(sample_map$accuratelocation))
mito_levels2 <- sort(unique(sample_map$mitotype))

hits <- vector("list", length(loc_levels) * length(mito_levels2))
ii <- 0L

for (loc in loc_levels) {
  loc_df <- sample_map %>% filter(accuratelocation == loc)

  for (m in mito_levels2) {
    s <- loc_df %>% filter(mitotype == m) %>% pull(sample)
    if (length(s) < 2) next  # can't define within-polymorphism with <2 samples

    poly_flag <- is_poly_within(allele_loc, s)
    idx <- which(poly_flag)
    if (length(idx) == 0) next

    ii <- ii + 1L
    hits[[ii]] <- tibble(
      variant_index = idx,
      accuratelocation = loc,
      mitotype = m
    )
  }
}

within_loc_mito_long <- bind_rows(hits[seq_len(ii)]) %>%
  left_join(var_anno, by = "variant_index") %>%
  filter(!is.na(gene), effect %in% c("syn", "nonsyn"))

within_loc_mito_pnps <- within_loc_mito_long %>%
  count(accuratelocation, mitotype, gene, effect, name = "n") %>%
  tidyr::complete(accuratelocation, mitotype, gene, effect, fill = list(n = 0L)) %>%
  tidyr::pivot_wider(names_from = effect, values_from = n, values_fill = 0L) %>%
  mutate(
    syn    = as.integer(syn),
    nonsyn = as.integer(nonsyn),
    total  = syn + nonsyn,
    pnps   = (nonsyn + PSEUDO_LOC) / (syn + PSEUDO_LOC),
    log10_pnps = log10(pnps)
  ) %>%
  filter(total > 0) %>%
  arrange(accuratelocation, gene, mitotype)

write_tsv(within_loc_mito_pnps, OUT_LOC_WITHIN_TSV)
log_msg("Wrote: ", OUT_LOC_WITHIN_TSV)

# Plot: genes on x, pn/ps on y, color by mitotype, facet by location
p_loc_within <- ggplot(within_loc_mito_pnps,
                       aes(x = gene, y = pnps, color = mitotype)) +
  geom_point(position = position_dodge(width = 0.6), size = 2.2, alpha = 0.9) +
  facet_wrap(~ accuratelocation, scales = "free_y") +
  scale_y_log10() +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "Gene",
    y = paste0("PN/PS (log scale; pseudocount ", PSEUDO_LOC, ")"),
    color = "Mitotype",
    title = "Within-location, within-mitotype PN/PS by gene (biallelic)",
    subtitle = "Within = polymorphic sites within each mitotype restricted to that location"
  )

ggsave(OUT_LOC_WITHIN_PNG, p_loc_within, width = 18, height = 10, dpi = 300)
log_msg("Wrote: ", OUT_LOC_WITHIN_PNG)





















# ============================================================
# (B) PAIRWISE BETWEEN-MITOTYPES: fixed differences proportions by gene
# variant contributes to pair if both mitotypes fixed and differ
# ============================================================
pair_df <- as.data.frame(t(combn(mito_levels, 2)), stringsAsFactors = FALSE)
colnames(pair_df) <- c("mitoA", "mitoB")

pair_hits <- vector("list", nrow(pair_df))
for (k in seq_len(nrow(pair_df))) {
  a <- pair_df$mitoA[k]
  b <- pair_df$mitoB[k]
  fa <- fixed_by_mito[, a]
  fb <- fixed_by_mito[, b]
  idx <- which(!is.na(fa) & !is.na(fb) & (fa != fb))
  if (length(idx) == 0) next

  pair_hits[[k]] <- var_anno[idx, , drop = FALSE] %>%
    mutate(mitoA = a, mitoB = b, pair = paste0(a, " vs ", b))
}

pair_long <- bind_rows(pair_hits) %>% filter(!is.na(gene))

pair_props <- pair_long %>%
  count(pair, mitoA, mitoB, gene, effect, name="n") %>%
  pivot_wider(names_from = effect, values_from = n, values_fill = 0) %>%
  mutate(
    total = syn + nonsyn,
    prop_nonsyn = ifelse(total > 0, nonsyn / total, NA_real_),
    prop_syn    = ifelse(total > 0, syn    / total, NA_real_)
  ) %>%
  filter(total >= MIN_TOTAL_PAIR) %>%
  arrange(gene, pair)

write_tsv(pair_props, OUT_PAIR_TSV)
log_msg("Wrote: ", OUT_PAIR_TSV)

p_pair <- ggplot(pair_props, aes(x = prop_nonsyn, y = prop_syn, color = pair)) +
  geom_point(aes(size = total), alpha = 0.8) +
  coord_equal(xlim = c(0,1), ylim = c(0,1)) +
  facet_wrap(~ gene) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  labs(
    x = "Proportion nonsynonymous (between pair)",
    y = "Proportion synonymous (between pair)",
    size = "N (fixed diffs)",
    color = "Pair",
    title = "Pairwise between-mitotype syn vs nonsyn proportions by gene (biallelic)"
  )

ggsave(OUT_PAIR_PNG, p_pair, width = 16, height = 11, dpi = 300)
log_msg("Wrote: ", OUT_PAIR_PNG)

# ============================================================
# (C) WITHIN-MITOTYPE: stacked counts (syn vs nonsyn) by gene
# ============================================================
within_counts <- within_long %>%
  count(mitotype, gene, effect, name="n") %>%
  complete(mitotype, gene, effect, fill=list(n=0)) %>%
  arrange(mitotype, gene, effect)

write_tsv(within_counts, OUT_BAR_TSV)
log_msg("Wrote: ", OUT_BAR_TSV)

p_bar <- ggplot(within_counts, aes(x = gene, y = n, fill = effect)) +
  geom_col() +
  facet_wrap(~ mitotype, scales = "free_y") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "Gene",
    y = "Within-mitotype polymorphic sites (count)",
    fill = "Effect",
    title = "Within-mitotype polymorphism counts: syn vs nonsyn (biallelic)"
  )

ggsave(OUT_BAR_PNG, p_bar, width = 14, height = 8, dpi = 300)
log_msg("Wrote: ", OUT_BAR_PNG)

log_msg("DONE")










# ----------------------------
# PN/PS settings
# ----------------------------
PSEUDO <- 1  # pseudocount to avoid Inf when PS=0

OUT_WITHIN_PNPS_TSV <- file.path(OUT_DIR, paste0(PREFIX, ".within.PNPS.tsv"))
OUT_WITHIN_PNPS_PNG <- file.path(OUT_DIR, paste0(PREFIX, ".within.PNPS.png"))

OUT_PAIR_PNPS_TSV   <- file.path(OUT_DIR, paste0(PREFIX, ".pairwise.between.PNPS.tsv"))
OUT_PAIR_PNPS_PNG   <- file.path(OUT_DIR, paste0(PREFIX, ".pairwise.between.PNPS.png"))

# ============================================================
# (D) WITHIN-MITOTYPE: PN/PS by gene
# ============================================================
within_pnps <- within_props %>%
  mutate(
    pnps = (nonsyn + PSEUDO) / (syn + PSEUDO),
    log10_pnps = log10(pnps)
  ) %>%
  arrange(gene, mitotype)

readr::write_tsv(within_pnps, OUT_WITHIN_PNPS_TSV)

p_within_pnps <- ggplot(within_pnps, aes(x = gene, y = pnps, color = mitotype)) +
  geom_point(position = position_dodge(width = 0.5), size = 2, alpha = 0.9) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "Gene",
    y = paste0("PN/PS ( (nonsyn+", PSEUDO, ")/(syn+", PSEUDO, ") )"),
    color = "Mitotype",
    title = "Within-mitotype PN/PS by gene (biallelic)"
  )

ggsave(OUT_WITHIN_PNPS_PNG, p_within_pnps, width = 14, height = 6, dpi = 300)

# Optional: log-scale version (often nicer)
OUT_WITHIN_LOGPNPS_PNG <- file.path(OUT_DIR, paste0(PREFIX, ".within.log10PNPS.png"))
p_within_logpnps <- ggplot(within_pnps, aes(x = gene, y = log10_pnps, color = mitotype)) +
  geom_point(position = position_dodge(width = 0.5), size = 2, alpha = 0.9) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "Gene",
    y = paste0("log10 PN/PS  (pseudocount ", PSEUDO, ")"),
    color = "Mitotype",
    title = "Within-mitotype log10(PN/PS) by gene (biallelic)"
  )
ggsave(OUT_WITHIN_LOGPNPS_PNG, p_within_logpnps, width = 14, height = 6, dpi = 300)



if (!exists("mito_levels")) mito_levels <- colnames(fixed_by_mito)
stopifnot(all(mito_levels %in% colnames(fixed_by_mito)))

# Variant annotation table aligned to fixed_by_mito rows
var_anno <- fix2 %>%
  transmute(
    variant_index = dplyr::row_number(),
    gene   = gene,
    effect = effect
  ) %>%
  filter(effect %in% c("syn","nonsyn"), !is.na(gene))

# All unordered pairs of mitotypes
pair_df <- as.data.frame(t(combn(mito_levels, 2)), stringsAsFactors = FALSE)
colnames(pair_df) <- c("mitoA", "mitoB")

pair_counts <- vector("list", nrow(pair_df))

for (k in seq_len(nrow(pair_df))) {
  a <- pair_df$mitoA[k]
  b <- pair_df$mitoB[k]

  fa <- fixed_by_mito[, a]
  fb <- fixed_by_mito[, b]

  # Variant is a between-difference for this pair if both fixed and different
  keep <- !is.na(fa) & !is.na(fb) & (fa != fb)
  idx  <- which(keep)
  if (length(idx) == 0) next

  tmp <- var_anno[idx, , drop = FALSE] %>%
    mutate(mitoA = a, mitoB = b)

  pair_counts[[k]] <- tmp
}

pair_long <- bind_rows(pair_counts) %>%
  mutate(pair = paste0(mitoA, " vs ", mitoB)) %>%
  filter(effect %in% c("syn","nonsyn"))

pair_props <- pair_long %>%
  count(pair, mitoA, mitoB, gene, effect, name = "n") %>%
  tidyr::pivot_wider(names_from = effect, values_from = n, values_fill = 0) %>%
  mutate(
    syn    = ifelse(is.na(syn), 0L, syn),
    nonsyn = ifelse(is.na(nonsyn), 0L, nonsyn),
    total  = syn + nonsyn,
    prop_nonsyn = ifelse(total > 0, nonsyn / total, NA_real_),
    prop_syn    = ifelse(total > 0, syn    / total, NA_real_)
  ) %>%
  filter(total > 0) %>%
  arrange(gene, pair)

# Save pairwise proportions table
readr::write_tsv(pair_props, OUT_PAIR_PNPS_TSV)

# ============================================================
# PN/PS plots for PAIRWISE BETWEEN
# ============================================================

pair_pnps <- pair_props %>%
  mutate(
    pnps = (nonsyn + PSEUDO) / (syn + PSEUDO),
    log10_pnps = log10(pnps)
  )

OUT_PAIR_PNPS_TSV   <- file.path(OUT_DIR, paste0(PREFIX, ".pairwise.between.PNPS.tsv"))
OUT_PAIR_PNPS_PNG   <- file.path(OUT_DIR, paste0(PREFIX, ".pairwise.between.PNPS.png"))
OUT_PAIR_LOGPNPS_PNG <- file.path(OUT_DIR, paste0(PREFIX, ".pairwise.between.log10PNPS.png"))

readr::write_tsv(pair_pnps, OUT_PAIR_PNPS_TSV)

p_pair_pnps <- ggplot(pair_pnps, aes(x = gene, y = pnps, color = pair, size = total)) +
  geom_point(alpha = 0.85, position = position_dodge(width = 0.6)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "Gene",
    y = paste0("PN/PS ( (nonsyn+", PSEUDO, ")/(syn+", PSEUDO, ") )"),
    color = "Pair",
    size = "N (fixed diffs)",
    title = "Pairwise between-mitotype PN/PS by gene (biallelic)"
  )

ggsave(OUT_PAIR_PNPS_PNG, p_pair_pnps, width = 16, height = 7, dpi = 300)

p_pair_logpnps <- ggplot(pair_pnps, aes(x = gene, y = log10_pnps, color = pair, size = total)) +
  geom_point(alpha = 0.85, position = position_dodge(width = 0.6)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "Gene",
    y = paste0("log10 PN/PS (pseudocount ", PSEUDO, ")"),
    color = "Pair",
    size = "N (fixed diffs)",
    title = "Pairwise between-mitotype log10(PN/PS) by gene (biallelic)"
  )

ggsave(OUT_PAIR_LOGPNPS_PNG, p_pair_logpnps, width = 16, height = 7, dpi = 300)

message("Wrote:\n", OUT_PAIR_PNPS_TSV, "\n", OUT_PAIR_PNPS_PNG, "\n", OUT_PAIR_LOGPNPS_PNG)




OUT_PAIR_PNPS_SCATTER_PNG <- file.path(
  OUT_DIR, paste0(PREFIX, ".pairwise.between.pN_vs_pS.png")
)

pair_pN_pS <- pair_props %>%
  mutate(
    pN = (nonsyn + PSEUDO) / (total + 2 * PSEUDO),
    pS = (syn    + PSEUDO) / (total + 2 * PSEUDO)
  )

# Sanity check (optional)
stopifnot(all(abs(pair_pN_pS$pN + pair_pN_pS$pS - 1) < 1e-6 | is.na(pair_pN_pS$pN)))

p_pN_pS <- ggplot(pair_pN_pS,
                  aes(x = pN, y = pS,
                      color = pair,
                      size = total)) +
  geom_point(alpha = 0.85) +
  facet_wrap(~ gene) +
  coord_equal(xlim = c(0,1), ylim = c(0,1)) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "pN = nonsyn / (syn + nonsyn)",
    y = "pS = syn / (syn + nonsyn)",
    color = "Mitotype pair",
    size  = "N (fixed differences)",
    title = "Pairwise between-mitotype pN vs pS by gene",
    subtitle = "Each point = one mitotype pair; biallelic sites only"
  )

ggsave(
  OUT_PAIR_PNPS_SCATTER_PNG,
  p_pN_pS,
  width = 14,
  height = 10,
  dpi = 300
)

message("Wrote pN vs pS plot:\n", OUT_PAIR_PNPS_SCATTER_PNG)

















# ----------------------------
# Mitotype subsets
# ----------------------------
SUBSET1 <- c("A","B","C","F","D")
SUBSET2 <- c("A","D","E","F","G","H")

# Keep within-mitotype rows
filter_within_set <- function(df, mitos) {
  df %>% dplyr::filter(mitotype %in% mitos)
}

# Keep pairwise rows: both endpoints inside set
filter_pair_set <- function(df, mitos) {
  df %>% dplyr::filter(mitoA %in% mitos, mitoB %in% mitos)
}


plot_within_counts <- function(df_counts, out_png, title_suffix="") {
  p <- ggplot(df_counts, aes(x = gene, y = n, fill = effect)) +
    geom_col() +
    facet_wrap(~ mitotype, scales = "free_y") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank()
    ) +
    labs(
      x = "Gene",
      y = "Within-mitotype polymorphic sites (count)",
      fill = "Effect",
      title = paste0("Within-mitotype syn vs nonsyn (by gene)", title_suffix)
    )

  ggsave(out_png, p, width = 14, height = 8, dpi = 300)
}

within_counts_s1 <- filter_within_set(within_counts, SUBSET1)
within_counts_s2 <- filter_within_set(within_counts, SUBSET2)

plot_within_counts(within_counts_s1,
  file.path(OUT_DIR, paste0(PREFIX, ".within.counts.ABCFD.png")),
  " — subset: A,B,C,F,D"
)

plot_within_counts(within_counts_s2,
  file.path(OUT_DIR, paste0(PREFIX, ".within.counts.ADEFGH.png")),
  " — subset: A,D,E,F,G,H"
)


plot_within_pN_pS <- function(df_props, out_png, title_suffix="") {
  p <- ggplot(df_props, aes(x = prop_nonsyn, y = prop_syn, color = mitotype)) +
    geom_point(size = 2, alpha = 0.9) +
    coord_equal(xlim = c(0,1), ylim = c(0,1)) +
    facet_wrap(~ gene) +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    labs(
      x = "pN = nonsyn/(syn+nonsyn)",
      y = "pS = syn/(syn+nonsyn)",
      color = "Mitotype",
      title = paste0("Within-mitotype pN vs pS by gene", title_suffix)
    )

  ggsave(out_png, p, width = 14, height = 10, dpi = 300)
}

within_props_s1 <- filter_within_set(within_props, SUBSET1)
within_props_s2 <- filter_within_set(within_props, SUBSET2)

plot_within_pN_pS(within_props_s1,
  file.path(OUT_DIR, paste0(PREFIX, ".within.pN_vs_pS.ABCFD.png")),
  " — subset: A,B,C,F,D"
)

plot_within_pN_pS(within_props_s2,
  file.path(OUT_DIR, paste0(PREFIX, ".within.pN_vs_pS.ADEFGH.png")),
  " — subset: A,D,E,F,G,H"
)


plot_within_pnps <- function(df, out_png, ycol="pnps", ylab="", title_suffix="") {
  p <- ggplot(df, aes(x = gene, y = .data[[ycol]], color = mitotype)) +
    geom_point(position = position_dodge(width = 0.5), size = 2, alpha = 0.9) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank()
    ) +
    labs(
      x = "Gene",
      y = ylab,
      color = "Mitotype",
      title = paste0("Within-mitotype ", ylab, " by gene", title_suffix)
    )

  ggsave(out_png, p, width = 14, height = 6, dpi = 300)
}

within_pnps_s1 <- filter_within_set(within_pnps, SUBSET1)
within_pnps_s2 <- filter_within_set(within_pnps, SUBSET2)

plot_within_pnps(within_pnps_s1,
  file.path(OUT_DIR, paste0(PREFIX, ".within.PNPS.ABCFD.png")),
  ycol="pnps",
  ylab=paste0("PN/PS (pseudocount ", PSEUDO, ")"),
  title_suffix=" — subset: A,B,C,F,D"
)

plot_within_pnps(within_pnps_s2,
  file.path(OUT_DIR, paste0(PREFIX, ".within.PNPS.ADEFGH.png")),
  ycol="pnps",
  ylab=paste0("PN/PS (pseudocount ", PSEUDO, ")"),
  title_suffix=" — subset: A,D,E,F,G,H"
)

plot_within_pnps(within_pnps_s1,
  file.path(OUT_DIR, paste0(PREFIX, ".within.log10PNPS.ABCFD.png")),
  ycol="log10_pnps",
  ylab=paste0("log10(PN/PS) (pseudocount ", PSEUDO, ")"),
  title_suffix=" — subset: A,B,C,F,D"
)

plot_within_pnps(within_pnps_s2,
  file.path(OUT_DIR, paste0(PREFIX, ".within.log10PNPS.ADEFGH.png")),
  ycol="log10_pnps",
  ylab=paste0("log10(PN/PS) (pseudocount ", PSEUDO, ")"),
  title_suffix=" — subset: A,D,E,F,G,H"
)





















GROUP1 <- c("A","B","C","F")  # merge these into one group

OUT_WITHIN_TSV <- file.path(OUT_DIR, paste0(PREFIX, ".within.PNPS.tsv"))
OUT_BETWEEN_TSV <- file.path(OUT_DIR, paste0(PREFIX, ".between.PNPS.tsv"))
OUT_COMBINED_TSV <- file.path(OUT_DIR, paste0(PREFIX, ".within_between.PNPS.tsv"))

OUT_WITHIN_PLOT <- file.path(OUT_DIR, paste0(PREFIX, ".within.PNPS.png"))
OUT_BETWEEN_PLOT <- file.path(OUT_DIR, paste0(PREFIX, ".between.PNPS.png"))
OUT_COMPARE_PLOT <- file.path(OUT_DIR, paste0(PREFIX, ".within_vs_between.PNPS.png"))

# ============================================================
# (1) Load mitotype map: expect 2 columns sample, mitotype (or force it)
# ============================================================
mito_df_raw <- suppressMessages(readr::read_delim(
  MITO_TSV,
  delim = ifelse(grepl("\\.csv$", MITO_TSV, ignore.case = TRUE), ",", "\t"),
  show_col_types = FALSE
))

# If file already has sample/mitotype, keep; else if it's 2-col unnamed, force names.
if (!all(c("sample","mitotype") %in% names(mito_df_raw))) {
  if (ncol(mito_df_raw) >= 2) {
    mito_df_raw <- mito_df_raw[, 1:2]
    colnames(mito_df_raw) <- c("sample","mitotype")
  } else {
    stop("MITO_TSV must have at least 2 columns: sample, mitotype")
  }
}

mito_df <- mito_df_raw %>%
  transmute(sample = as.character(sample),
            mitotype = as.character(mitotype)) %>%
  filter(!is.na(sample), !is.na(mitotype)) %>%
  distinct()

# Merge mitotypes into two groups
mito_df <- mito_df %>%
  mutate(group = ifelse(mitotype %in% GROUP1, "ABCF", "Other"))

# ============================================================
# (2) Read VCF + keep BIALLELIC only (post-split or not)
# ============================================================
vcf <- vcfR::read.vcfR(VCF_IN, verbose = FALSE)
gt_all <- vcfR::extract.gt(vcf, element = "GT")  # variants x samples
fix_all <- as.data.frame(vcf@fix, stringsAsFactors = FALSE)

# Keep biallelic only: ALT has no comma, and not missing
bial_mask <- !is.na(fix_all$ALT) & !str_detect(fix_all$ALT, ",")
fix <- fix_all[bial_mask, , drop = FALSE]
gt  <- gt_all[bial_mask, , drop = FALSE]

message("Biallelic sites kept: ", nrow(fix))

# Keep only samples present in mito map
keep_samples <- intersect(colnames(gt), mito_df$sample)
if (length(keep_samples) < 2) stop("Too few VCF samples match mitotype map. Check sample naming.")
gt <- gt[, keep_samples, drop = FALSE]
mito_df <- mito_df %>% filter(sample %in% keep_samples)

# Group sample indices by merged group
group_levels <- c("ABCF", "Other")
group_to_samples <- split(mito_df$sample, mito_df$group)
group_to_samples <- group_to_samples[group_levels]  # keep order if both exist

if (any(vapply(group_to_samples, length, integer(1)) < 2)) {
  message("WARNING: One group has <2 samples; within-polymorphism for that group may be sparse or undefined.")
}

# ============================================================
# (3) Parse ANN -> gene + syn/nonsyn (robust vapply)
# ============================================================
parse_ann <- function(info) {
  if (is.na(info) || !str_detect(info, "ANN=")) return(c(NA_character_, NA_character_))

  ann <- str_replace(info, "^.*ANN=", "")
  ann <- str_split(ann, ";", simplify = TRUE)[1]
  entries <- str_split(ann, ",")[[1]]

  # Prefer nonsyn if any entry has missense; else syn if any has synonymous
  hit_gene_nonsyn <- NA_character_
  hit_gene_syn <- NA_character_

  for (e in entries) {
    f <- str_split(e, "\\|", simplify = TRUE)
    if (ncol(f) < 5) next
    eff  <- f[1,2]
    gene <- f[1,4]
    if (is.na(gene) || gene == "") next

    if (str_detect(eff, "missense_variant")) {
      hit_gene_nonsyn <- gene
      break
    }
    if (is.na(hit_gene_syn) && str_detect(eff, "synonymous_variant")) {
      hit_gene_syn <- gene
    }
  }

  if (!is.na(hit_gene_nonsyn)) return(c(hit_gene_nonsyn, "nonsyn"))
  if (!is.na(hit_gene_syn))    return(c(hit_gene_syn, "syn"))
  c(NA_character_, NA_character_)
}

info_vec <- vcfR::getINFO(vcf)[bial_mask]  # INFO for biallelic rows
anno_mat <- t(vapply(info_vec, parse_ann, FUN.VALUE = character(2)))
fix$gene   <- anno_mat[,1]
fix$effect <- anno_mat[,2]

# Keep only syn/nonsyn with gene
fix2 <- fix %>% filter(effect %in% c("syn","nonsyn"), !is.na(gene))
if (nrow(fix2) == 0) stop("No syn/nonsyn annotations found in biallelic sites. Check ANN parsing / VCF.")

# Subset GT to those rows (need row indices in the biallelic-filtered matrix)
# fix2 rownames are original rownames from fix_all subset; they should match row positions in gt after masking.
# Because we subset via bial_mask, rownames(fix) are preserved from fix_all; gt rows align with fix by position.
# So use row positions in fix (not original). We'll map using rownames in fix.
row_pos_in_fix <- match(rownames(fix2), rownames(fix))
stopifnot(all(!is.na(row_pos_in_fix)))

gt2 <- gt[row_pos_in_fix, , drop = FALSE]
stopifnot(nrow(gt2) == nrow(fix2))

# Convert GT to allele index (haploid or diploid; take first allele)
gt_to_allele <- function(x) {
  if (is.na(x) || x %in% c(".", "./.", ".|.")) return(NA_integer_)
  a <- str_split(x, "[/|]", simplify = TRUE)[1]
  suppressWarnings(as.integer(a))
}

first_allele <- function(x) {
  x <- sub("[/|].*$", "", x)          # keep first allele
  x[x %in% c(".", "./.", ".|.")] <- NA
  suppressWarnings(as.integer(x))
}

allele_mat <- matrix(first_allele(gt2),
                     nrow = nrow(gt2),
                     ncol = ncol(gt2),
                     dimnames = dimnames(gt2))
# ============================================================
# (4) Classify variants as within-polymorphic per group; and between fixed differences
# ============================================================
fixed_by_group <- matrix(NA_integer_, nrow=nrow(allele_mat), ncol=length(group_levels),
                         dimnames=list(NULL, group_levels))
poly_by_group  <- matrix(FALSE, nrow=nrow(allele_mat), ncol=length(group_levels),
                         dimnames=list(NULL, group_levels))

for (j in seq_along(group_levels)) {
  g <- group_levels[j]
  s <- group_to_samples[[g]]
  if (is.null(s) || length(s) == 0) next
  sub <- allele_mat[, s, drop = FALSE]

  aset <- apply(sub, 1, function(v) sort(unique(v[!is.na(v)])))
  poly <- vapply(aset, function(a) length(a) > 1, logical(1))
  fixed <- vapply(aset, function(a) if (length(a)==1) a[1] else NA_integer_, integer(1))

  poly_by_group[, j]  <- poly
  fixed_by_group[, j] <- fixed
}

# Between-group fixed difference: both groups fixed and different
fa <- fixed_by_group[, "ABCF"]
fb <- fixed_by_group[, "Other"]
between_flags <- !is.na(fa) & !is.na(fb) & (fa != fb)

# ============================================================
# (5) Summarise syn/nonsyn counts and PN/PS
# ============================================================
var_anno <- fix2 %>% transmute(variant_index = seq_len(nrow(fix2)), gene, effect)

# Within: variants that are polymorphic within each group
within_df <- as.data.frame(poly_by_group) %>%
  mutate(variant_index = seq_len(nrow(poly_by_group))) %>%
  pivot_longer(-variant_index, names_to="group", values_to="is_within") %>%
  filter(is_within) %>%
  select(variant_index, group)

within_anno <- within_df %>%
  left_join(var_anno, by="variant_index") %>%
  filter(!is.na(gene), effect %in% c("syn","nonsyn"))

within_counts <- within_anno %>%
  count(group, gene, effect, name="n") %>%
  complete(group, gene, effect, fill = list(n=0L)) %>%
  pivot_wider(names_from = effect, values_from = n, values_fill = 0L) %>%
  mutate(
    syn = as.integer(syn),
    nonsyn = as.integer(nonsyn),
    total = syn + nonsyn,
    pnps = (nonsyn + PSEUDO) / (syn + PSEUDO),
    log10_pnps = log10(pnps),
    category = "Within (polymorphic)"
  ) %>%
  arrange(gene, group)

# Between: fixed differences between the two groups (one row per variant)
between_idx <- which(between_flags)
between_anno <- var_anno[between_idx, , drop = FALSE] %>%
  mutate(groupA="ABCF", groupB="Other") %>%
  filter(!is.na(gene), effect %in% c("syn","nonsyn"))

between_counts <- between_anno %>%
  count(gene, effect, name="n") %>%
  complete(gene, effect, fill = list(n=0L)) %>%
  pivot_wider(names_from = effect, values_from = n, values_fill = 0L) %>%
  mutate(
    syn = as.integer(syn),
    nonsyn = as.integer(nonsyn),
    total = syn + nonsyn,
    pnps = (nonsyn + PSEUDO) / (syn + PSEUDO),
    log10_pnps = log10(pnps),
    category = "Between (fixed diffs)",
    group = "ABCF vs Other"
  ) %>%
  arrange(gene)

# Write TSVs
readr::write_tsv(within_counts, OUT_WITHIN_TSV)
readr::write_tsv(between_counts, OUT_BETWEEN_TSV)

combined_pnps <- bind_rows(
  within_counts %>% select(category, group, gene, syn, nonsyn, total, pnps, log10_pnps),
  between_counts %>% select(category, group, gene, syn, nonsyn, total, pnps, log10_pnps)
)
readr::write_tsv(combined_pnps, OUT_COMBINED_TSV)

# ============================================================
# (6) Plots
# ============================================================

# Within PN/PS by gene (points; colored by group)
p_within <- ggplot(within_counts %>% filter(total > 0),
                   aes(x = gene, y = pnps, color = group)) +
  geom_point(position = position_dodge(width = 0.5), size = 2.2, alpha = 0.9) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank()) +
  labs(
    x = "Gene",
    y = paste0("PN/PS = (nonsyn+", PSEUDO, ")/(syn+", PSEUDO, ")"),
    color = "Group",
    title = "Within-group PN/PS by gene (ABCF merged vs Other)",
    subtitle = "Within = polymorphic variants within each group (biallelic only)"
  )
ggsave(OUT_WITHIN_PLOT, p_within, width = 14, height = 6, dpi = 300)

# Between PN/PS by gene (single group)
p_between <- ggplot(between_counts %>% filter(total > 0),
                    aes(x = gene, y = pnps)) +
  geom_point(size = 2.4, alpha = 0.9) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank()) +
  labs(
    x = "Gene",
    y = paste0("PN/PS = (nonsyn+", PSEUDO, ")/(syn+", PSEUDO, ")"),
    title = "Between-group PN/PS by gene (ABCF vs Other)",
    subtitle = "Between = fixed differences between groups (biallelic only)"
  )
ggsave(OUT_BETWEEN_PLOT, p_between, width = 14, height = 6, dpi = 300)

# Compare within vs between on the same axis (log10 PN/PS)
p_compare <- ggplot(combined_pnps %>% filter(total > 0),
                    aes(x = gene, y = log10_pnps, color = category, shape = group)) +
  geom_point(position = position_dodge(width = 0.6), size = 2.2, alpha = 0.9) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank()) +
  labs(
    x = "Gene",
    y = paste0("log10(PN/PS) (pseudocount ", PSEUDO, ")"),
    color = "Category",
    shape = "Group",
    title = "Within vs Between log10(PN/PS) by gene (ABCF merged vs Other)",
    subtitle = "Within = polymorphic; Between = fixed differences; biallelic only"
  )
ggsave(OUT_COMPARE_PLOT, p_compare, width = 16, height = 7, dpi = 300)

message("Wrote:\n",
        OUT_WITHIN_TSV, "\n",
        OUT_BETWEEN_TSV, "\n",
        OUT_COMBINED_TSV, "\n",
        OUT_WITHIN_PLOT, "\n",
        OUT_BETWEEN_PLOT, "\n",
        OUT_COMPARE_PLOT)










OUT_PAIR_BETWEEN_TSV <- file.path(OUT_DIR, paste0(PREFIX, ".between_fixed.PNPS.tsv"))
OUT_PAIR_BETWEEN_PNG <- file.path(OUT_DIR, paste0(PREFIX, ".between_fixed.PNPS.png"))

OUT_PAIR_WITHIN_TSV  <- file.path(OUT_DIR, paste0(PREFIX, ".within_poly.PNPS.tsv"))
OUT_PAIR_WITHIN_PNG  <- file.path(OUT_DIR, paste0(PREFIX, ".within_poly.PNPS.png"))




group_levels <- c("ABCF","Other")
group_to_samples <- split(mito_df$sample, mito_df$group)[group_levels]

fixed_by_group <- matrix(NA_integer_, nrow=nrow(allele_mat), ncol=2, dimnames=list(NULL, group_levels))
for (j in seq_along(group_levels)) {
  g <- group_levels[j]
  s <- group_to_samples[[g]]
  sub <- allele_mat[, s, drop = FALSE]
  aset <- apply(sub, 1, function(v) sort(unique(v[!is.na(v)])))
  fixed_by_group[, j] <- vapply(aset, function(a) if (length(a)==1) a[1] else NA_integer_, integer(1))
}

fa <- fixed_by_group[, "ABCF"]
fb <- fixed_by_group[, "Other"]
between_idx <- which(!is.na(fa) & !is.na(fb) & (fa != fb))

pair_between <- var_anno[between_idx, , drop = FALSE] %>%
  mutate(pair = "ABCF vs Other") %>%
  count(pair, gene, effect, name="n") %>%
  complete(pair, gene, effect, fill=list(n=0L)) %>%
  pivot_wider(names_from = effect, values_from = n, values_fill = 0L) %>%
  mutate(
    syn = as.integer(syn),
    nonsyn = as.integer(nonsyn),
    total = syn + nonsyn,
    pN = nonsyn,
    pS = syn,
    pnps = (nonsyn + PSEUDO) / (syn + PSEUDO),
    log10_pnps = log10(pnps)
  ) %>%
  filter(total > 0) %>%
  arrange(gene)

readr::write_tsv(pair_between, OUT_PAIR_BETWEEN_TSV)

p_between <- ggplot(pair_between, aes(x = pN, y = pS,col= gene)) +
  geom_point(size = 2.6, alpha = 0.9) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  labs(
    x = "pN (count of nonsyn fixed differences)",
    y = "pS (count of syn fixed differences)",
    title = "Between merged groups (ABCF vs Other): pN vs pS by gene (biallelic)",
    subtitle = "Each facet is a gene; one point per pair (only one pair here)"
  )

ggsave(OUT_PAIR_BETWEEN_PNG, p_between, width = 14, height = 10, dpi = 300)

# ============================================================
# (4B) Pairwise WITHIN (polymorphism) across pairs of ORIGINAL mitotypes,
# then re-labeled by merged-group pair: ABCF vs ABCF, Other vs Other, ABCF vs Other.
# For each original pair, count variants where BOTH mitotypes are fixed but DIFFER (a divergence-like within pair),
# OR (optional) count variants polymorphic within either mitotype (more "within" style).
#
# Below: "within-pair divergence": both mitotypes fixed and different.
# This produces MANY points per gene.
# ============================================================
pair_mat <- t(combn(mito_levels, 2))
pair_df <- tibble(mitoA = pair_mat[,1], mitoB = pair_mat[,2]) %>%
  left_join(mito_df %>% distinct(mitotype, group), by = c("mitoA" = "mitotype")) %>%
  rename(groupA = group) %>%
  left_join(mito_df %>% distinct(mitotype, group), by = c("mitoB" = "mitotype")) %>%
  rename(groupB = group) %>%
  mutate(
    merged_pair = ifelse(groupA == groupB, paste0(groupA, " vs ", groupB),
                         "ABCF vs Other")
  )

# fixed allele per ORIGINAL mitotype
fixed_by_mito <- matrix(NA_integer_, nrow=nrow(allele_mat), ncol=length(mito_levels),
                        dimnames=list(NULL, mito_levels))
for (m in mito_levels) {
  s <- mito_to_samples[[m]]
  sub <- allele_mat[, s, drop = FALSE]
  aset <- apply(sub, 1, function(v) sort(unique(v[!is.na(v)])))
  fixed_by_mito[, m] <- vapply(aset, function(a) if (length(a)==1) a[1] else NA_integer_, integer(1))
}

pair_counts <- vector("list", nrow(pair_df))

for (k in seq_len(nrow(pair_df))) {
  a <- pair_df$mitoA[k]
  b <- pair_df$mitoB[k]

  fa <- fixed_by_mito[, a]
  fb <- fixed_by_mito[, b]

  idx <- which(!is.na(fa) & !is.na(fb) & (fa != fb))
  if (length(idx) == 0) next

  pair_counts[[k]] <- var_anno[idx, , drop = FALSE] %>%
    mutate(
      mitoA = a, mitoB = b,
      pair = paste0(a, " vs ", b),
      merged_pair = pair_df$merged_pair[k]
    )
}

pair_long <- bind_rows(pair_counts)
if (nrow(pair_long) == 0) stop("No pairwise fixed-difference sites found among original mitotypes (biallelic).")

pair_within <- pair_long %>%
  count(merged_pair, pair, mitoA, mitoB, gene, effect, name="n") %>%
  tidyr::pivot_wider(names_from = effect, values_from = n, values_fill = 0L) %>%
  mutate(
    syn = ifelse(is.na(syn), 0L, as.integer(syn)),
    nonsyn = ifelse(is.na(nonsyn), 0L, as.integer(nonsyn)),
    total = syn + nonsyn,
    pN = nonsyn,
    pS = syn,
    pnps = (nonsyn + PSEUDO) / (syn + PSEUDO),
    log10_pnps = log10(pnps)
  ) %>%
  filter(total > 0) %>%
  arrange(gene, merged_pair, pair)

readr::write_tsv(pair_within, OUT_PAIR_WITHIN_TSV)

# Plot: for each gene, each ORIGINAL pair is a point; colored by merged_pair group label
p_within_pairs <- ggplot(pair_within, aes(x = pN, y = pS, color = merged_pair)) +
  geom_point(aes(size = total), alpha = 0.85) +
  facet_wrap(~ gene, scales = "free") +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  labs(
    x = "pN (count of nonsyn differences)",
    y = "pS (count of syn differences)",
    color = "Merged pair label",
    size = "N (pN+pS)",
    title = "Pairwise comparisons among original mitotypes: pN vs pS by gene (biallelic)",
    subtitle = "Each point = one mitotype pair; points colored by merged group-pair (ABCF vs Other, ABCF vs ABCF, Other vs Other)"
  )

ggsave(OUT_PAIR_WITHIN_PNG, p_within_pairs, width = 16, height = 10, dpi = 300)

message("Wrote:\n",
        OUT_PAIR_BETWEEN_TSV, "\n", OUT_PAIR_BETWEEN_PNG, "\n",
        OUT_PAIR_WITHIN_TSV, "\n", OUT_PAIR_WITHIN_PNG)











OUT_PAIR_BETWEEN_TSV2 <- file.path(OUT_DIR, paste0(PREFIX, ".between_fixed2.PNPS.tsv"))
OUT_PAIR_BETWEEN_PNG2 <- file.path(OUT_DIR, paste0(PREFIX, ".between_fixed2.PNPS.png"))


pair_between <- var_anno[between_idx, , drop = FALSE] %>%
  mutate(pair = "ABCF vs Other") %>%
  count(pair, gene, effect, name = "n") %>%
  complete(pair, gene, effect, fill = list(n = 0L)) %>%
  pivot_wider(names_from = effect, values_from = n, values_fill = 0L) %>%
  mutate(
    syn    = as.integer(syn),
    nonsyn = as.integer(nonsyn),
    total  = syn + nonsyn,

    # ---- percentages ----
    pN_pct = nonsyn / total,
    pS_pct = syn    / total,

    # optional PN/PS still retained
    pnps = (nonsyn + PSEUDO) / (syn + PSEUDO),
    log10_pnps = log10(pnps)
  ) %>%
  filter(total > 0) %>%
  arrange(gene)

readr::write_tsv(pair_between, OUT_PAIR_BETWEEN_TSV2)

# -------------------------
# Plot: percent pN vs percent pS
# -------------------------
p_between <- ggplot(pair_between, aes(x = pN_pct, y = pS_pct, color = gene)) +
  geom_point(size = 2.6, alpha = 0.9) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  labs(
    x = "pN (%) – nonsynonymous fixed differences",
    y = "pS (%) – synonymous fixed differences",
    color = "Gene",
    title = "Between merged groups (ABCF vs Other): pN% vs pS% by gene (biallelic)",
    subtitle = "Percentages within each gene; pN% + pS% = 100"
  )

ggsave(OUT_PAIR_BETWEEN_PNG2, p_between, width = 14, height = 10, dpi = 300)


















GROUP1 <- c("A","B","C","F")
mito_df2 <- mito_df %>%
  mutate(group = ifelse(mitotype %in% GROUP1, "ABCF", "Other"))

# Ensure sample order aligns with allele_mat columns
stopifnot(all(colnames(allele_mat) %in% mito_df2$sample))
mito_df2 <- mito_df2 %>% filter(sample %in% colnames(allele_mat))

samps_abcf  <- mito_df2 %>% filter(group == "ABCF")  %>% pull(sample)
samps_other <- mito_df2 %>% filter(group == "Other") %>% pull(sample)

stopifnot(length(samps_abcf)  >= 2)
stopifnot(length(samps_other) >= 2)

# Helper: for a sample set, mark variants polymorphic within set
is_poly_within <- function(allele_mat, samples) {
  sub <- allele_mat[, samples, drop = FALSE]
  apply(sub, 1, function(v) {
    vv <- v[!is.na(v)]
    length(unique(vv)) > 1
  })
}

# Helper: fixed allele within set (NA if not fixed or no calls)
fixed_allele <- function(allele_mat, samples) {
  sub <- allele_mat[, samples, drop = FALSE]
  apply(sub, 1, function(v) {
    vv <- v[!is.na(v)]
    u <- unique(vv)
    if (length(u) == 1) u[1] else NA_integer_
  })
}

# =========================
# (1) WITHIN ABCF
# =========================
within_abcf_idx <- which(is_poly_within(allele_mat, samps_abcf))
within_abcf <- var_anno[within_abcf_idx, , drop = FALSE] %>%
  mutate(comparison = "Within ABCF") %>%
  count(comparison, gene, effect, name="n") %>%
  complete(comparison, gene, effect, fill=list(n=0L)) %>%
  pivot_wider(names_from = effect, values_from = n, values_fill = 0L) %>%
  mutate(
    syn    = as.integer(syn),
    nonsyn = as.integer(nonsyn),
    total  = syn + nonsyn,
    pN_pct = 100 * nonsyn / total,
    pS_pct = 100 * syn    / total,
    pnps   = (nonsyn + PSEUDO) / (syn + PSEUDO),
    log10_pnps = log10(pnps)
  ) %>%
  filter(total > 0) %>%
  arrange(gene)

OUT_WITHIN_ABCF_TSV <- file.path(OUT_DIR, paste0(PREFIX, ".within.ABCF.pN_pS_pct.tsv"))
OUT_WITHIN_ABCF_PNG <- file.path(OUT_DIR, paste0(PREFIX, ".within.ABCF.pN_pS_pct.png"))
readr::write_tsv(within_abcf, OUT_WITHIN_ABCF_TSV)

p_within_abcf <- ggplot(within_abcf, aes(x = pN_pct, y = pS_pct, color = gene)) +
  geom_point(size = 2.6, alpha = 0.9) +
  coord_equal(xlim = c(0,100), ylim = c(0,100)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  labs(
    x = "pN (%) – nonsynonymous (within ABCF polymorphisms)",
    y = "pS (%) – synonymous (within ABCF polymorphisms)",
    color = "Gene",
    title = "Within ABCF: pN% vs pS% by gene (biallelic)",
    subtitle = "Only sites polymorphic within ABCF; pN% + pS% = 100"
  )
ggsave(OUT_WITHIN_ABCF_PNG, p_within_abcf, width = 12, height = 7, dpi = 300)

# =========================
# (2) WITHIN Other
# =========================
within_other_idx <- which(is_poly_within(allele_mat, samps_other))
within_other <- var_anno[within_other_idx, , drop = FALSE] %>%
  mutate(comparison = "Within Other") %>%
  count(comparison, gene, effect, name="n") %>%
  complete(comparison, gene, effect, fill=list(n=0L)) %>%
  pivot_wider(names_from = effect, values_from = n, values_fill = 0L) %>%
  mutate(
    syn    = as.integer(syn),
    nonsyn = as.integer(nonsyn),
    total  = syn + nonsyn,
    pN_pct = 100 * nonsyn / total,
    pS_pct = 100 * syn    / total,
    pnps   = (nonsyn + PSEUDO) / (syn + PSEUDO),
    log10_pnps = log10(pnps)
  ) %>%
  filter(total > 0) %>%
  arrange(gene)

OUT_WITHIN_OTHER_TSV <- file.path(OUT_DIR, paste0(PREFIX, ".within.Other.pN_pS_pct.tsv"))
OUT_WITHIN_OTHER_PNG <- file.path(OUT_DIR, paste0(PREFIX, ".within.Other.pN_pS_pct.png"))
readr::write_tsv(within_other, OUT_WITHIN_OTHER_TSV)

p_within_other <- ggplot(within_other, aes(x = pN_pct, y = pS_pct, color = gene)) +
  geom_point(size = 2.6, alpha = 0.9) +
  coord_equal(xlim = c(0,100), ylim = c(0,100)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  labs(
    x = "pN (%) – nonsynonymous (within Other polymorphisms)",
    y = "pS (%) – synonymous (within Other polymorphisms)",
    color = "Gene",
    title = "Within Other: pN% vs pS% by gene (biallelic)",
    subtitle = "Only sites polymorphic within Other; pN% + pS% = 100"
  )
ggsave(OUT_WITHIN_OTHER_PNG, p_within_other, width = 12, height = 7, dpi = 300)

# =========================
# (3) BETWEEN ABCF vs Other (fixed differences)
# =========================
fa <- fixed_allele(allele_mat, samps_abcf)
fb <- fixed_allele(allele_mat, samps_other)

between_idx <- which(!is.na(fa) & !is.na(fb) & (fa != fb))

between_abcf_other <- var_anno[between_idx, , drop = FALSE] %>%
  mutate(comparison = "ABCF vs Other (between)") %>%
  count(comparison, gene, effect, name="n") %>%
  complete(comparison, gene, effect, fill=list(n=0L)) %>%
  pivot_wider(names_from = effect, values_from = n, values_fill = 0L) %>%
  mutate(
    syn    = as.integer(syn),
    nonsyn = as.integer(nonsyn),
    total  = syn + nonsyn,
    pN_pct = 100 * nonsyn / total,
    pS_pct = 100 * syn    / total,
    pnps   = (nonsyn + PSEUDO) / (syn + PSEUDO),
    log10_pnps = log10(pnps)
  ) %>%
  filter(total > 0) %>%
  arrange(gene)

OUT_BETWEEN_TSV <- file.path(OUT_DIR, paste0(PREFIX, ".between.ABCF_vs_Other.pN_pS_pct.tsv"))
OUT_BETWEEN_PNG <- file.path(OUT_DIR, paste0(PREFIX, ".between.ABCF_vs_Other.pN_pS_pct.png"))
readr::write_tsv(between_abcf_other, OUT_BETWEEN_TSV)

p_between <- ggplot(between_abcf_other, aes(x = pN_pct, y = pS_pct, color = gene)) +
  geom_point(size = 2.6, alpha = 0.9) +
  coord_equal(xlim = c(0,100), ylim = c(0,100)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  labs(
    x = "pN (%) – nonsynonymous fixed differences",
    y = "pS (%) – synonymous fixed differences",
    color = "Gene",
    title = "Between ABCF vs Other: pN% vs pS% by gene (biallelic)",
    subtitle = "Only sites fixed in ABCF and fixed in Other, with different fixed alleles; pN% + pS% = 100"
  )
ggsave(OUT_BETWEEN_PNG, p_between, width = 12, height = 7, dpi = 300)

message("Wrote:\n",
        OUT_WITHIN_ABCF_TSV, "\n", OUT_WITHIN_ABCF_PNG, "\n",
        OUT_WITHIN_OTHER_TSV, "\n", OUT_WITHIN_OTHER_PNG, "\n",
        OUT_BETWEEN_TSV, "\n", OUT_BETWEEN_PNG)


















samps_abcf  <- mito_df2 %>% filter(group == "ABCF")  %>% pull(sample)
samps_other <- mito_df2 %>% filter(group == "Other") %>% pull(sample)
stopifnot(length(samps_abcf) >= 2, length(samps_other) >= 2)

# ----------------------------
# Helpers
# ----------------------------
is_poly_within <- function(allele_mat, samples) {
  sub <- allele_mat[, samples, drop = FALSE]
  apply(sub, 1, function(v) {
    vv <- v[!is.na(v)]
    length(unique(vv)) > 1
  })
}

fixed_allele <- function(allele_mat, samples) {
  sub <- allele_mat[, samples, drop = FALSE]
  apply(sub, 1, function(v) {
    vv <- v[!is.na(v)]
    u <- unique(vv)
    if (length(u) == 1) u[1] else NA_integer_
  })
}

summarise_gene_props <- function(df, label) {
  df %>%
    mutate(comparison = label) %>%
    count(comparison, gene, effect, name = "n") %>%
    complete(comparison, gene, effect, fill = list(n = 0L)) %>%
    pivot_wider(names_from = effect, values_from = n, values_fill = 0L) %>%
    mutate(
      syn    = as.integer(syn),
      nonsyn = as.integer(nonsyn),
      total  = syn + nonsyn,
      pN = ifelse(total > 0, nonsyn / total, NA_real_),
      pS = ifelse(total > 0, syn    / total, NA_real_),
      pnps = (nonsyn + PSEUDO) / (syn + PSEUDO),
      log10_pnps = log10(pnps)
    ) %>%
    filter(total > 0) %>%
    arrange(gene)
}

# ============================================================
# (1) WITHIN ABCF: pN vs pS (0..1)
# ============================================================
within_abcf_idx <- which(is_poly_within(allele_mat, samps_abcf))
within_abcf_tbl <- summarise_gene_props(var_anno[within_abcf_idx, , drop = FALSE], "Within ABCF")

OUT_WITHIN_ABCF_TSV <- file.path(OUT_DIR, paste0(PREFIX, ".within.ABCF.pN_pS.tsv"))
OUT_WITHIN_ABCF_PNG <- file.path(OUT_DIR, paste0(PREFIX, ".within.ABCF.pN_pS.png"))
write_tsv(within_abcf_tbl, OUT_WITHIN_ABCF_TSV)

p_within_abcf <- ggplot(within_abcf_tbl, aes(x = pN, y = pS, color = gene)) +
  geom_point(size = 2.6, alpha = 0.9) +
  coord_equal(xlim = c(0,1), ylim = c(0,1)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  labs(
    x = "pN (nonsyn / (syn+nonsyn))",
    y = "pS (syn / (syn+nonsyn))",
    color = "Gene",
    title = "Within ABCF: pN vs pS by gene (biallelic)"
  )
ggsave(OUT_WITHIN_ABCF_PNG, p_within_abcf, width = 12, height = 7, dpi = 300)

# ============================================================
# (2) WITHIN Other: pN vs pS (0..1)
# ============================================================
within_other_idx <- which(is_poly_within(allele_mat, samps_other))
within_other_tbl <- summarise_gene_props(var_anno[within_other_idx, , drop = FALSE], "Within Other")

OUT_WITHIN_OTHER_TSV <- file.path(OUT_DIR, paste0(PREFIX, ".within.Other.pN_pS.tsv"))
OUT_WITHIN_OTHER_PNG <- file.path(OUT_DIR, paste0(PREFIX, ".within.Other.pN_pS.png"))
write_tsv(within_other_tbl, OUT_WITHIN_OTHER_TSV)

p_within_other <- ggplot(within_other_tbl, aes(x = pN, y = pS, color = gene)) +
  geom_point(size = 2.6, alpha = 0.9) +
  coord_equal(xlim = c(0,1), ylim = c(0,1)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  labs(
    x = "pN (nonsyn / (syn+nonsyn))",
    y = "pS (syn / (syn+nonsyn))",
    color = "Gene",
    title = "Within Other: pN vs pS by gene (biallelic)"
  )
ggsave(OUT_WITHIN_OTHER_PNG, p_within_other, width = 12, height = 7, dpi = 300)

# ============================================================
# (3) BETWEEN ABCF vs Other: fixed diffs pN vs pS (0..1)
# ============================================================
fa <- fixed_allele(allele_mat, samps_abcf)
fb <- fixed_allele(allele_mat, samps_other)
between_idx <- which(!is.na(fa) & !is.na(fb) & (fa != fb))

between_tbl <- summarise_gene_props(var_anno[between_idx, , drop = FALSE], "ABCF vs Other (between)")

OUT_BETWEEN_TSV <- file.path(OUT_DIR, paste0(PREFIX, ".between.ABCF_vs_Other.pN_pS.tsv"))
OUT_BETWEEN_PNG <- file.path(OUT_DIR, paste0(PREFIX, ".between.ABCF_vs_Other.pN_pS.png"))
write_tsv(between_tbl, OUT_BETWEEN_TSV)

p_between <- ggplot(between_tbl, aes(x = pN, y = pS, color = gene)) +
  geom_point(size = 2.6, alpha = 0.9) +
  coord_equal(xlim = c(0,1), ylim = c(0,1)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  labs(
    x = "pN (nonsyn / (syn+nonsyn))",
    y = "pS (syn / (syn+nonsyn))",
    color = "Gene",
    title = "Between ABCF vs Other: pN vs pS by gene (biallelic)"
  )
ggsave(OUT_BETWEEN_PNG, p_between, width = 12, height = 7, dpi = 300)

# ============================================================
# (4) NEW: Genes on x, PN/PS by mitotype (A,B,C,...) within each GROUP
#     Here "within mitotype" means polymorphic within that mitotype.
# ============================================================

# per-variant fixed/polymorphic status within each mitotype
mito_levels <- sort(unique(mito_df2$mitotype))
mito_to_samples <- split(mito_df2$sample, mito_df2$mitotype)

poly_by_mito <- matrix(FALSE, nrow = nrow(allele_mat), ncol = length(mito_levels),
                       dimnames = list(NULL, mito_levels))

for (m in mito_levels) {
  s <- mito_to_samples[[m]]
  if (length(s) < 2) next
  poly_by_mito[, m] <- is_poly_within(allele_mat, s)
}

within_mito_long <- as.data.frame(poly_by_mito) %>%
  mutate(variant_index = seq_len(nrow(poly_by_mito))) %>%
  pivot_longer(-variant_index, names_to = "mitotype", values_to = "is_within") %>%
  filter(is_within) %>%
  select(variant_index, mitotype)

# attach gene/effect
var_anno2 <- var_anno %>%
  mutate(variant_index = seq_len(nrow(var_anno)))

within_mito_anno <- within_mito_long %>%
  left_join(var_anno2, by = "variant_index") %>%
  filter(effect %in% c("syn","nonsyn"), !is.na(gene))

# add group labels
within_mito_anno <- within_mito_anno %>%
  left_join(mito_df2 %>% distinct(mitotype, group), by = "mitotype")

# summarise per (group, mitotype, gene)
within_mito_pnps <- within_mito_anno %>%
  count(group, mitotype, gene, effect, name = "n") %>%
  complete(group, mitotype, gene, effect, fill = list(n = 0L)) %>%
  pivot_wider(names_from = effect, values_from = n, values_fill = 0L) %>%
  mutate(
    syn    = as.integer(syn),
    nonsyn = as.integer(nonsyn),
    total  = syn + nonsyn,
    pnps   = (nonsyn + PSEUDO) / (syn + PSEUDO),
    log10_pnps = log10(pnps)
  ) %>%
  filter(total > 0) %>%
  arrange(group, mitotype, gene)

OUT_WITHIN_MITO_PNPS_TSV <- file.path(OUT_DIR, paste0(PREFIX, ".within_mitotype.PNPS_by_gene.tsv"))
OUT_WITHIN_MITO_PNPS_PNG <- file.path(OUT_DIR, paste0(PREFIX, ".within_mitotype.PNPS_by_gene.png"))
write_tsv(within_mito_pnps, OUT_WITHIN_MITO_PNPS_TSV)

p_mito_pnps <- ggplot(within_mito_pnps, aes(x = gene, y = pnps, color = mitotype)) +
  geom_point(position = position_dodge(width = 0.6), alpha = 0.9) +
  facet_wrap(~ group, scales = "free_y") +
  theme_bw() + ylim(0, 10)+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "Gene",
    y = paste0("PN/PS ( (nonsyn+", PSEUDO, ")/(syn+", PSEUDO, ") )"),
    color = "Mitotype",
    title = "Within-mitotype PN/PS by gene, faceted by merged group (biallelic)"
  )
ggsave(OUT_WITHIN_MITO_PNPS_PNG, p_mito_pnps, width = 16, height = 7, dpi = 300)

message("Wrote:\n",
        OUT_WITHIN_ABCF_TSV, "\n", OUT_WITHIN_ABCF_PNG, "\n",
        OUT_WITHIN_OTHER_TSV, "\n", OUT_WITHIN_OTHER_PNG, "\n",
        OUT_BETWEEN_TSV, "\n", OUT_BETWEEN_PNG, "\n",
        OUT_WITHIN_MITO_PNPS_TSV, "\n", OUT_WITHIN_MITO_PNPS_PNG)


















samples <- read.csv("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/Sample_ID_Species_merged_20251227.csv")

usobtusasamps <- subset(samples, Species == "Daphnia obtusa")


usobtusasamps$preferred_label <- ifelse(
  !is.na(usobtusasamps$Sample_ID_old) & usobtusasamps$Sample_ID_old != "",
  usobtusasamps$Sample_ID_old,
  usobtusasamps$Sample_ID
)






metadata2 <- metadata %>%
  mutate(clone = str_remove(clone, "_clone"))


## ------------------------
## Tip table: label + mitotype
## Label rule: Sample_ID_old > accuratelocation > clone
## ------------------------

usobtusasamps_1 <- usobtusasamps %>%
  distinct(preferred_label, .keep_all = TRUE)

samps <- tibble(clone = usobtusasamps_1$preferred_label) %>%
  left_join(metadata2, by = "clone") %>%
  left_join(usobtusasamps_1, by = c("clone" = "preferred_label"))

samps <- as.data.frame(samps)
samps <- samps %>%
  mutate(
    accuratelocation = if_else(
      str_detect(clone, "^(EBG|FS|PYR|bdw|RAP)"),
      str_extract(clone, "^(EBG|FS|PYR|bdw|RAP)"),
      accuratelocation
    )
  )



samps2 <- samps %>%
  select(clone, accuratelocation, mitotype) %>%
  mutate(
    accuratelocation = na_if(accuratelocation, ""),
    mitotype = if_else(is.na(mitotype) | mitotype == "", "Unknown", mitotype)
  )

# Optional: drop Unknown mitotypes or missing locations
samps2 <- samps2 %>% filter(!is.na(accuratelocation), mitotype != "Unknown")












PSEUDO_LOC=1

OUT_LOC_PNPS_TSV <- file.path(OUT_DIR, paste0(PREFIX, ".location_mitotype.PNPS.tsv"))
OUT_LOC_PNPS_PNG <- file.path(OUT_DIR, paste0(PREFIX, ".location_mitotype.PNPS.png"))

# ---- 1) sample -> mitotype (A/B/C/...) using mitotypes table (CloneA, Group)
# If mitotypes is already in memory:
stopifnot(exists("mitotypes"))
mito_key <- mitotypes %>%
  transmute(sample = as.character(CloneA),
            mitotype = as.character(Group)) %>%
  distinct()

# Merge with mito_df (VCF samples that were kept)
sample_info <- mito_df %>%
  transmute(sample = as.character(sample)) %>%
  left_join(mito_key, by = "sample")

# ---- 2) sample -> accuratelocation
# Prefer metadata if you have it (clone column often has _clone suffix)
if (exists("metadata")) {
  meta_loc <- metadata %>%
    transmute(
      clone_key = str_remove(as.character(clone), "_clone$"),
      accuratelocation = na_if(as.character(accuratelocation), "")
    ) %>%
    distinct()

  sample_info <- sample_info %>%
    mutate(clone_key = sample) %>%   # sample names are clone-like IDs
    left_join(meta_loc, by = "clone_key") %>%
    select(-clone_key)
} else {
  # fallback: prefix-based location labels
  sample_info <- sample_info %>%
    mutate(
      accuratelocation = if_else(
        str_detect(sample, regex("^(EBG|FS|PYR|bdw|RAP)", ignore_case = TRUE)),
        str_to_upper(str_extract(sample, regex("^(EBG|FS|PYR|bdw|RAP)", ignore_case = TRUE))),
        NA_character_
      )
    )
}

# Keep only samples with both mitotype and location
sample_info <- sample_info %>%
  mutate(
    mitotype = if_else(is.na(mitotype) | mitotype == "", "Unknown", mitotype),
    accuratelocation = na_if(accuratelocation, "")
  ) %>%
  filter(!is.na(accuratelocation), mitotype != "Unknown")

# Subset allele_mat columns to these samples (and keep ordering)
keep_samps_loc <- intersect(colnames(allele_mat), sample_info$sample)
allele_loc <- allele_mat[, keep_samps_loc, drop = FALSE]
sample_info <- sample_info %>% filter(sample %in% keep_samps_loc)

log_msg("Location+mitotype samples kept: ", nrow(sample_info),
        " across ", n_distinct(sample_info$accuratelocation), " locations and ",
        n_distinct(sample_info$mitotype), " mitotypes.")

# ---- helpers
is_poly_within_samples <- function(allele_mat, samples) {
  if (length(samples) < 2) return(rep(FALSE, nrow(allele_mat)))
  sub <- allele_mat[, samples, drop = FALSE]
  apply(sub, 1, function(v) {
    vv <- v[!is.na(v)]
    length(unique(vv)) > 1
  })
}

fixed_allele_samples <- function(allele_mat, samples) {
  if (length(samples) < 1) return(rep(NA_integer_, nrow(allele_mat)))
  sub <- allele_mat[, samples, drop = FALSE]
  apply(sub, 1, function(v) {
    vv <- v[!is.na(v)]
    u <- unique(vv)
    if (length(u) == 1) u[1] else NA_integer_
  })
}

# ---- 3) Compute within-polymorphic variants per (location, mitotype)
loc_levels  <- sort(unique(sample_info$accuratelocation))
mito_levels_loc <- sort(unique(sample_info$mitotype))

within_hits <- list()
between_hits <- list()

for (loc in loc_levels) {
  loc_samples <- sample_info %>% filter(accuratelocation == loc)

  for (m in mito_levels_loc) {
    s_m  <- loc_samples %>% filter(mitotype == m) %>% pull(sample)
    s_ot <- loc_samples %>% filter(mitotype != m) %>% pull(sample)

    # --- WITHIN: polymorphic within mitotype at this location
    poly_flag <- is_poly_within_samples(allele_loc, s_m)
    idx_within <- which(poly_flag)
    if (length(idx_within) > 0) {
      within_hits[[length(within_hits) + 1]] <-
        var_anno[idx_within, , drop = FALSE] %>%
        mutate(
          accuratelocation = loc,
          mitotype = m,
          comparison = "Within (poly within mitotype)"
        )
    }

    # --- BETWEEN: fixed differences m vs others at this location
    fa <- fixed_allele_samples(allele_loc, s_m)
    fb <- fixed_allele_samples(allele_loc, s_ot)
    idx_between <- which(!is.na(fa) & !is.na(fb) & (fa != fb))
    if (length(idx_between) > 0) {
      between_hits[[length(between_hits) + 1]] <-
        var_anno[idx_between, , drop = FALSE] %>%
        mutate(
          accuratelocation = loc,
          mitotype = m,
          comparison = "Between (fixed m vs other mitotypes)"
        )
    }
  }
}

loc_long <- bind_rows(within_hits, between_hits) %>%
  filter(!is.na(gene), effect %in% c("syn", "nonsyn"))

# ---- 4) Summarise pN/pS per location x mitotype x comparison (optionally per gene too)
loc_pnps <- loc_long %>%
  count(comparison, accuratelocation, mitotype, effect, name = "n") %>%
  complete(comparison, accuratelocation, mitotype, effect, fill = list(n = 0L)) %>%
  pivot_wider(names_from = effect, values_from = n, values_fill = 0L) %>%
  mutate(
    syn = as.integer(syn),
    nonsyn = as.integer(nonsyn),
    total = syn + nonsyn,
    pnps = (nonsyn + PSEUDO_LOC) / (syn + PSEUDO_LOC),
    log10_pnps = log10(pnps)
  ) %>%
  filter(total > 0)

write_tsv(loc_pnps, OUT_LOC_PNPS_TSV)
log_msg("Wrote: ", OUT_LOC_PNPS_TSV)

# ---- 5) Plot by location, color by mitotype, facet by comparison
p_loc <- ggplot(loc_pnps,
                aes(x = accuratelocation, y = pnps, color = mitotype)) +
  geom_point(position = position_dodge(width = 0.6), size = 2.2, alpha = 0.9) +
  scale_y_log10() +
  facet_wrap(~ comparison, ncol = 1, scales = "free_y") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "Location",
    y = paste0("PN/PS (log scale; pseudocount ", PSEUDO_LOC, ")"),
    color = "Mitotype",
    title = "PN/PS by location and mitotype",
    subtitle = "Within = polymorphic within mitotype at location; Between = fixed differences vs other mitotypes at location"
  )

ggsave(OUT_LOC_PNPS_PNG, p_loc, width = 14, height = 8, dpi = 300)
log_msg("Wrote: ", OUT_LOC_PNPS_PNG)
