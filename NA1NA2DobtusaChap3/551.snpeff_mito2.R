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


PSEUDO_LOC=0.1
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
meta_loc <- samps %>%
  transmute(
    sample = str_remove(as.character(Sample_ID), "_clone$"),
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


within_loc_mito_pnps <- within_loc_mito_pnps %>%
  mutate(label = paste0("pN=", nonsyn, ", pS=", syn))


write_tsv(within_loc_mito_pnps, OUT_LOC_WITHIN_TSV)
log_msg("Wrote: ", OUT_LOC_WITHIN_TSV)

# Plot: genes on x, pn/ps on y, color by mitotype, facet by location
p_loc_within <- ggplot(within_loc_mito_pnps,
                       aes(x = gene, y = pnps, color = mitotype)) +
  geom_point(position = position_dodge(width = 0.6), size = 2.2, alpha = 0.9) +
  geom_text(
    aes(label = label),
    position = position_dodge(width = 0.6),
    vjust = -0.6,
    size = 3,
    show.legend = FALSE
  ) +
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
    subtitle = "Point labels are pN:pS (nonsyn:syn)"
  )


ggsave(OUT_LOC_WITHIN_PNG, p_loc_within, width = 18, height = 10, dpi = 300)
log_msg("Wrote: ", OUT_LOC_WITHIN_PNG)


















# ============================================================
# (B) WITHIN EACH POPULATION (location): PN/PS by gene
#     "Within population" = sites polymorphic among ALL samples in that location
# ============================================================

OUT_POP_WITHIN_TSV <- file.path(OUT_DIR, paste0(PREFIX, ".within_population.PNPS_by_gene.tsv"))
OUT_POP_WITHIN_PNG <- file.path(OUT_DIR, paste0(PREFIX, ".within_population.PNPS_by_gene.png"))

# helper: polymorphic within a set of samples
is_poly_within <- function(allele_mat, samples) {
  if (length(samples) < 2) return(rep(FALSE, nrow(allele_mat)))
  sub <- allele_mat[, samples, drop = FALSE]
  apply(sub, 1, function(v) {
    vv <- v[!is.na(v)]
    length(unique(vv)) > 1
  })
}

pop_hits <- list()
ii <- 0L

for (loc in sort(unique(sample_map$accuratelocation))) {
  s_loc <- sample_map %>% filter(accuratelocation == loc) %>% pull(sample)
  if (length(s_loc) < 2) next

  idx <- which(is_poly_within(allele_loc, s_loc))
  if (length(idx) == 0) next

  ii <- ii + 1L
  pop_hits[[ii]] <- tibble(
    variant_index = idx,
    accuratelocation = loc
  )
}

within_pop_long <- bind_rows(pop_hits) %>%
  left_join(var_anno, by = "variant_index") %>%
  filter(!is.na(gene), effect %in% c("syn","nonsyn"))

within_pop_pnps <- within_pop_long %>%
  count(accuratelocation, gene, effect, name = "n") %>%
  complete(accuratelocation, gene, effect, fill = list(n = 0L)) %>%
  pivot_wider(names_from = effect, values_from = n, values_fill = 0L) %>%
  mutate(
    syn    = as.integer(syn),
    nonsyn = as.integer(nonsyn),
    total  = syn + nonsyn,
    pnps   = (nonsyn + PSEUDO_LOC) / (syn + PSEUDO_LOC),
    log10_pnps = log10(pnps),
    label = paste0("pN=", nonsyn, ", pS=", syn)
  ) %>%
  filter(total > 0) %>%
  arrange(accuratelocation, gene)

write_tsv(within_pop_pnps, OUT_POP_WITHIN_TSV)
log_msg("Wrote: ", OUT_POP_WITHIN_TSV)

p_within_pop <- ggplot(within_pop_pnps, aes(x = gene, y = pnps)) +
  geom_point(size = 2.2, alpha = 0.9) +
  geom_text(aes(label = label), vjust = -0.6, size = 3, show.legend = FALSE) +
  facet_wrap(~ accuratelocation, scales = "free_y") +
  scale_y_log10() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank()) +
  labs(
    x = "Gene",
    y = paste0("PN/PS (log scale; pseudocount ", PSEUDO_LOC, ")"),
    title = "Within-population PN/PS by gene (biallelic)",
    subtitle = "Within = polymorphic sites among all samples in that location"
  )

ggsave(OUT_POP_WITHIN_PNG, p_within_pop, width = 18, height = 10, dpi = 300)
log_msg("Wrote: ", OUT_POP_WITHIN_PNG)


# ============================================================
# (C) BETWEEN MITOTYPES WITHIN EACH POPULATION: PN/PS by gene
#     "Between" = fixed differences between mitotypes at that location
# ============================================================

OUT_POP_BETWEEN_TSV <- file.path(OUT_DIR, paste0(PREFIX, ".between_mitotypes_within_population.PNPS_by_gene.tsv"))
OUT_POP_BETWEEN_PNG <- file.path(OUT_DIR, paste0(PREFIX, ".between_mitotypes_within_population.PNPS_by_gene.png"))

# helper: fixed allele within a set of samples (NA if not fixed or no calls)
fixed_allele <- function(allele_mat, samples) {
  if (length(samples) < 1) return(rep(NA_integer_, nrow(allele_mat)))
  sub <- allele_mat[, samples, drop = FALSE]
  apply(sub, 1, function(v) {
    vv <- v[!is.na(v)]
    u <- unique(vv)
    if (length(u) == 1) u[1] else NA_integer_
  })
}

between_hits <- list()
kk <- 0L

for (loc in sort(unique(sample_map$accuratelocation))) {
  loc_df <- sample_map %>% filter(accuratelocation == loc)

  mitos_here <- sort(unique(loc_df$mitotype))
  if (length(mitos_here) < 2) next

  pairs <- combn(mitos_here, 2, simplify = FALSE)

  for (pr in pairs) {
    a <- pr[1]; b <- pr[2]
    s_a <- loc_df %>% filter(mitotype == a) %>% pull(sample)
    s_b <- loc_df %>% filter(mitotype == b) %>% pull(sample)
    if (length(s_a) < 1 || length(s_b) < 1) next

    fa <- fixed_allele(allele_loc, s_a)
    fb <- fixed_allele(allele_loc, s_b)

    idx <- which(!is.na(fa) & !is.na(fb) & fa != fb)
    if (length(idx) == 0) next

    kk <- kk + 1L
    between_hits[[kk]] <- tibble(
      variant_index = idx,
      accuratelocation = loc,
      mitoA = a,
      mitoB = b,
      pair = paste0(a, " vs ", b)
    )
  }
}

between_pop_long <- bind_rows(between_hits) %>%
  left_join(var_anno, by = "variant_index") %>%
  filter(!is.na(gene), effect %in% c("syn","nonsyn"))

between_pop_pnps <- between_pop_long %>%
  count(accuratelocation, pair, mitoA, mitoB, gene, effect, name = "n") %>%
  complete(accuratelocation, pair, mitoA, mitoB, gene, effect, fill = list(n = 0L)) %>%
  pivot_wider(names_from = effect, values_from = n, values_fill = 0L) %>%
  mutate(
    syn    = as.integer(syn),
    nonsyn = as.integer(nonsyn),
    total  = syn + nonsyn,
    pnps   = (nonsyn + PSEUDO_LOC) / (syn + PSEUDO_LOC),
    log10_pnps = log10(pnps),
    label = paste0("pN=", nonsyn, ", pS=", syn)
  ) %>%
  filter(total > 0) %>%
  arrange(accuratelocation, gene, pair)

write_tsv(between_pop_pnps, OUT_POP_BETWEEN_TSV)
log_msg("Wrote: ", OUT_POP_BETWEEN_TSV)

p_between_pop <- ggplot(between_pop_pnps,
                        aes(x = gene, y = pnps, color = pair)) +
  geom_point(position = position_dodge(width = 0.7), size = 2.2, alpha = 0.9) +
  geom_text(aes(label = label),
            position = position_dodge(width = 0.7),
            vjust = -0.6, size = 3, show.legend = FALSE) +
  facet_wrap(~ accuratelocation, scales = "free_y") +
  scale_y_log10() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank()) +
  labs(
    x = "Gene",
    y = paste0("PN/PS (log scale; pseudocount ", PSEUDO_LOC, ")"),
    color = "Mitotype pair",
    title = "Between-mitotype PN/PS by gene within each population (biallelic)",
    subtitle = "Between = fixed differences where both mitotypes are fixed (within that location)"
  )

ggsave(OUT_POP_BETWEEN_PNG, p_between_pop, width = 30, height = 10, dpi = 300)
log_msg("Wrote: ", OUT_POP_BETWEEN_PNG)

















plot_ps_vs_pn <- function(df, out_png, title, subtitle = NULL) {
  gg <- ggplot(df, aes(x = nonsyn, y = syn, color = pair)) +
    geom_point(size = 2.2, alpha = 0.9) +
    facet_grid(accuratelocation ~ gene, scales = "free") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 0),
      panel.grid.minor = element_blank(),
      strip.text.x = element_text(size = 8),
      strip.text.y = element_text(size = 8)
    ) +
    labs(
      x = "pN (nonsyn count)",
      y = "pS (syn count)",
      color = "Pair",
      title = title,
      subtitle = subtitle
    )

  ggsave(out_png, gg, width = 18, height = 10, dpi = 300)
  log_msg("Wrote: ", out_png)
  gg
}


within_mito_scatter <- within_loc_mito_pnps %>%
  transmute(
    accuratelocation,
    gene,
    pair = paste0("Within ", mitotype),
    syn,
    nonsyn
  ) %>%
  filter(syn + nonsyn > 0)

OUT_WITHIN_MITO_PS_PN <- file.path(OUT_DIR, paste0(PREFIX, ".within_mitotype.ps_vs_pn.by_gene_by_pop.png"))

p1 <- plot_ps_vs_pn(
  within_mito_scatter,
  OUT_WITHIN_MITO_PS_PN,
  title = "Within-mitotype (per population): pS vs pN by gene",
  subtitle = "Each point = one (population, gene, mitotype); counts from polymorphic sites within that mitotype"
)


within_pop_scatter <- within_pop_pnps %>%
  transmute(
    accuratelocation,
    gene,
    pair = "Within population",
    syn,
    nonsyn
  ) %>%
  filter(syn + nonsyn > 0)

OUT_WITHIN_POP_PS_PN <- file.path(OUT_DIR, paste0(PREFIX, ".within_population.ps_vs_pn.by_gene_by_pop.png"))

p2 <- plot_ps_vs_pn(
  within_pop_scatter,
  OUT_WITHIN_POP_PS_PN,
  title = "Within-population: pS vs pN by gene",
  subtitle = "Each point = one (population, gene); counts from polymorphic sites within the population"
)

between_pair_scatter <- between_pop_pnps %>%
  transmute(
    accuratelocation,
    gene,
    pair = pair,
    syn,
    nonsyn
  ) %>%
  filter(syn + nonsyn > 0)

OUT_BETWEEN_PAIR_PS_PN <- file.path(OUT_DIR, paste0(PREFIX, ".between_mitotype_pairs.ps_vs_pn.by_gene_by_pop.png"))

p3 <- plot_ps_vs_pn(
  between_pair_scatter,
  OUT_BETWEEN_PAIR_PS_PN,
  title = "Between mitotype pairs (within population): pS vs pN by gene",
  subtitle = "Each point = one (population, gene, mitotype-pair); counts from fixed differences between the two mitotypes"
)





OUT_ALL_PS_PN_TSV <- file.path(OUT_DIR, paste0(PREFIX, ".all_comparisons.ps_vs_pn.tsv"))
OUT_ALL_PS_PN_PNG <- file.path(OUT_DIR, paste0(PREFIX, ".all_comparisons.ps_vs_pn.by_gene.png"))

# ============================================================
# Combine all comparison types into one table
# ============================================================

df_pspn_all <- bind_rows(
  # (1) within mitotype (within pond)
  within_loc_mito_pnps %>%
    transmute(
      gene,
      accuratelocation,
      comparison = "Within mitotype",
      syn,
      nonsyn
    ),

  # (2) within pond (ignoring mitotype)
  within_pop_pnps %>%
    transmute(
      gene,
      accuratelocation,
      comparison = "Within population",
      syn,
      nonsyn
    ),

  # (3) between mitotype pairs within pond
  between_pop_pnps %>%
    transmute(
      gene,
      accuratelocation,
      comparison = "Between mitotypes (pairwise)",
      syn,
      nonsyn
    )
) %>%
  filter(!is.na(gene), !is.na(accuratelocation)) %>%
  mutate(
    comparison = factor(
      comparison,
      levels = c(
        "Within population",
        "Within mitotype",
        "Between mitotypes (pairwise)"
      )
    ),
    pond = factor(accuratelocation)
  ) %>%
  filter((syn + nonsyn) > 0)

write_tsv(df_pspn_all, OUT_ALL_PS_PN_TSV)
log_msg("Wrote: ", OUT_ALL_PS_PN_TSV)


p_all_pspn <- ggplot(
  df_pspn_all,
  aes(x = nonsyn, y = syn, color = pond, shape = comparison)
) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    linewidth = 0.6,
    color = "grey40"
  ) +
  geom_point(
    size = 2.6,
    alpha = 0.9,
    position = position_jitter(width = 0.05, height = 0.05)
  ) +
  facet_wrap(~ gene, scales = "free") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0),
    legend.position = "right"
  ) +
  labs(
    x = "pN (nonsyn count)",
    y = "pS (syn count)",
    color = "Pond",
    shape = "Comparison",
    title = "pS vs pN by gene",
    subtitle = "Dashed line shows pN = pS (slope = 1)"
  )


ggsave(OUT_ALL_PS_PN_PNG, p_all_pspn, width = 16, height = 10, dpi = 300)
log_msg("Wrote: ", OUT_ALL_PS_PN_PNG)



















# ============================================================
# Make 3 separate figures (same style):
# facet by gene, color by pond, slope=1 line
# ============================================================

plot_pspn_by_gene <- function(df, out_png, title, subtitle = NULL) {
  df <- df %>%
    filter(!is.na(gene), !is.na(accuratelocation)) %>%
    mutate(pond = factor(accuratelocation)) %>%
    filter((syn + nonsyn) > 0)

  p <- ggplot(df, aes(x = nonsyn, y = syn, color = pond)) +
    geom_abline(
      slope = 1, intercept = 0,
      linetype = "dashed", linewidth = 0.6, color = "grey40"
    ) +
    geom_point(
      size = 2.6, alpha = 0.9,
      position = position_jitter(width = 0.05, height = 0.05)
    ) +
    facet_wrap(~ gene, scales = "free") +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 0),
      legend.position = "right"
    ) +
    labs(
      x = "pN (nonsyn count)",
      y = "pS (syn count)",
      color = "Pond",
      title = title,
      subtitle = subtitle
    )

  ggsave(out_png, p, width = 16, height = 10, dpi = 300)
  log_msg("Wrote: ", out_png)
  p
}

# ----------------------------
# (1) Within mitotype (within pond)
# ----------------------------
df_within_mito <- within_loc_mito_pnps %>%
  transmute(gene, accuratelocation, syn, nonsyn)

OUT_WITHIN_MITO_PS_PN_PNG <- file.path(
  OUT_DIR, paste0(PREFIX, ".within_mitotype.ps_vs_pn.by_gene.png")
)

p_within_mito_pspn <- plot_pspn_by_gene(
  df_within_mito,
  OUT_WITHIN_MITO_PS_PN_PNG,
  title = "Within mitotype (per pond): pS vs pN by gene",
  subtitle = "Dashed line shows pN = pS (slope = 1)"
)

# ----------------------------
# (2) Within population (pond)
# ----------------------------
df_within_pop <- within_pop_pnps %>%
  transmute(gene, accuratelocation, syn, nonsyn)

OUT_WITHIN_POP_PS_PN_PNG <- file.path(
  OUT_DIR, paste0(PREFIX, ".within_population.ps_vs_pn.by_gene.png")
)

p_within_pop_pspn <- plot_pspn_by_gene(
  df_within_pop,
  OUT_WITHIN_POP_PS_PN_PNG,
  title = "Within pond: pS vs pN by gene",
  subtitle = "Dashed line shows pN = pS (slope = 1)"
)




df_between_pairs <- between_pop_pnps %>%
  transmute(
    gene,
    accuratelocation,
    pair,          # keep the mitotype comparison (e.g. "A vs B")
    syn,
    nonsyn
  ) %>%
  filter((syn + nonsyn) > 0) %>%
  mutate(
    pond = factor(accuratelocation),
    pair = factor(pair)
  )

OUT_BETWEEN_PAIRS_PS_PN_PNG <- file.path(
  OUT_DIR, paste0(PREFIX, ".between_mitotype_pairs.ps_vs_pn.by_gene.png")
)

p_between_pairs_pspn <- ggplot(
  df_between_pairs,
  aes(x = nonsyn, y = syn, color = pond, shape = pair)
) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    linewidth = 0.6,
    color = "grey40"
  ) +
  geom_point(
    size = 2.6,
    alpha = 0.9,
    position = position_jitter(width = 0.05, height = 0.05)
  ) +
  facet_wrap(~ gene, scales = "free") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right"
  ) + xlim(0,30) + ylim(0,80)
  labs(
    x = "pN (nonsyn count)",
    y = "pS (syn count)",
    color = "Pond",
    shape = "Mitotype pair",
    title = "Between mitotype pairs (within pond): pS vs pN by gene",
    subtitle = "Each point = one mitotype comparison within a pond; dashed line is pN = pS"
  )

ggsave(
  OUT_BETWEEN_PAIRS_PS_PN_PNG,
  p_between_pairs_pspn,
  width = 16,
  height = 10,
  dpi = 300
)

log_msg("Wrote: ", OUT_BETWEEN_PAIRS_PS_PN_PNG)





















PSEUDO_PAIR <- if (exists("PSEUDO_LOC")) PSEUDO_LOC else 0.1

OUT_PAIRWISE_ALL_TSV <- file.path(OUT_DIR, paste0(PREFIX, ".pairwise_all_mitotypes.PNPS_by_gene.tsv"))
OUT_PAIRWISE_ALL_PNG <- file.path(OUT_DIR, paste0(PREFIX, ".pairwise_all_mitotypes.ps_vs_pn.by_gene.png"))

# Sanity checks
stopifnot(exists("mito_levels"), exists("fixed_by_mito"), exists("var_anno"))
stopifnot(all(c("variant_index","gene","effect") %in% names(var_anno)))
stopifnot(nrow(fixed_by_mito) == nrow(var_anno))

# All unordered pairs of mitotypes
pair_mat <- t(combn(mito_levels, 2))
pair_df  <- tibble(mitoA = pair_mat[, 1], mitoB = pair_mat[, 2]) %>%
  mutate(pair = paste0(mitoA, " vs ", mitoB))

# Collect per-pair variant hits (fixed in both mitotypes and different)
pair_hits <- vector("list", nrow(pair_df))

for (k in seq_len(nrow(pair_df))) {
  a <- pair_df$mitoA[k]
  b <- pair_df$mitoB[k]

  fa <- fixed_by_mito[, a]
  fb <- fixed_by_mito[, b]

  idx <- which(!is.na(fa) & !is.na(fb) & (fa != fb))
  if (length(idx) == 0) next

  pair_hits[[k]] <- var_anno[idx, , drop = FALSE] %>%
    transmute(
      pair = pair_df$pair[k],
      mitoA = a,
      mitoB = b,
      gene,
      effect
    )
}

pair_long_all <- bind_rows(pair_hits) %>%
  filter(!is.na(gene), effect %in% c("syn","nonsyn"))

if (nrow(pair_long_all) == 0) stop("No fixed differences found between any mitotype pairs (with current fixed_by_mito definition).")

# Summarise syn/nonsyn counts per (pair, gene)
pairwise_all <- pair_long_all %>%
  count(pair, mitoA, mitoB, gene, effect, name = "n") %>%
  tidyr::complete(pair, mitoA, mitoB, gene, effect, fill = list(n = 0L)) %>%
  tidyr::pivot_wider(names_from = effect, values_from = n, values_fill = 0L) %>%
  mutate(
    syn    = as.integer(syn),
    nonsyn = as.integer(nonsyn),
    total  = syn + nonsyn,
    pN = nonsyn,
    pS = syn,
    pnps = (nonsyn + PSEUDO_PAIR) / (syn + PSEUDO_PAIR),
    log10_pnps = log10(pnps)
  ) %>%
  filter(total > 0) %>%
  arrange(gene, pair)

write_tsv(pairwise_all, OUT_PAIRWISE_ALL_TSV)
log_msg("Wrote: ", OUT_PAIRWISE_ALL_TSV)

# Plot: raw pS vs pN, colored by pair, faceted by gene, slope=1 line
p_pairwise_all <- ggplot(pairwise_all, aes(x = pN, y = pS, color = pair)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.6, color = "grey40") +
  geom_point(size = 2.4, alpha = 0.85, position = position_jitter(width = 0.05, height = 0.05)) +
  facet_wrap(~ gene, scales = "free") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0),
    legend.position = "right"
  ) +
  labs(
    x = "pN (nonsyn fixed-diff count)",
    y = "pS (syn fixed-diff count)",
    color = "Mitotype pair",
    title = "Pairwise between all mitotypes: pS vs pN by gene (ignoring population)",
    subtitle = "Each point = one mitotype pair; variants counted when both mitotypes are fixed and differ; dashed line is pN = pS"
  )

ggsave(OUT_PAIRWISE_ALL_PNG, p_pairwise_all, width = 16, height = 10, dpi = 300)
log_msg("Wrote: ", OUT_PAIRWISE_ALL_PNG)



NA2 <- c("A", "B", "C", "F")


pairwise_all <- pairwise_all %>%
  mutate(
    groupA = if_else(mitoA %in% NA2, "NA2", "NA1"),
    groupB = if_else(mitoB %in% NA2, "NA2", "NA1"),

    comparison_group = case_when(
      groupA == "NA2" & groupB == "NA2" ~ "NA2",
      groupA == "NA1" & groupB == "NA1" ~ "NA1",
      TRUE                              ~ "Between NA1 and NA2"
    ),

    comparison_group = factor(
      comparison_group,
      levels = c("NA2", "NA1", "Between NA1 and NA2")
    )
  )



OUT_PAIRWISE_ALL_PNG <- file.path(
  OUT_DIR, paste0(PREFIX, ".pairwise_all_mitotypes.ps_vs_pn.by_gene.NA_groups.png")
)

p_pairwise_all <- ggplot(
  pairwise_all,
  aes(x = pN, y = pS, color = pair, shape = comparison_group)
) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    linewidth = 0.6,
    color = "grey40"
  ) +
  geom_point(
    size = 2.6,
    alpha = 0.85,
    position = position_jitter(width = 0.05, height = 0.05)
  ) +
  facet_wrap(~ gene, scales = "free") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0),
    legend.position = "right"
  ) +
  labs(
    x = "pN (nonsyn fixed-difference count)",
    y = "pS (syn fixed-difference count)",
    color = "Mitotype pair",
    shape = "Comparison group",
    title = "Pairwise pS vs pN between all mitotypes (by gene)",
    subtitle = "Shapes indicate NA groups; dashed line shows pN = pS"
  )

ggsave(OUT_PAIRWISE_ALL_PNG, p_pairwise_all, width = 16, height = 10, dpi = 300)
log_msg("Wrote: ", OUT_PAIRWISE_ALL_PNG)





ggsave(
  sub("\\.png$", ".pdf", OUT_WITHIN_PNG),
  p_within,
  width = 14,
  height = 10,
  device = cairo_pdf
)

# 2) Within-location, within-mitotype PN/PS
ggsave(
  sub("\\.png$", ".pdf", OUT_LOC_WITHIN_PNG),
  p_loc_within,
  width = 18,
  height = 10,
  device = cairo_pdf
)

# 3) Within-population PN/PS
ggsave(
  sub("\\.png$", ".pdf", OUT_POP_WITHIN_PNG),
  p_within_pop,
  width = 18,
  height = 10,
  device = cairo_pdf
)

# 4) Between-mitotypes within population PN/PS
ggsave(
  sub("\\.png$", ".pdf", OUT_POP_BETWEEN_PNG),
  p_between_pop,
  width = 30,
  height = 10,
  device = cairo_pdf
)

# 5) Within-mitotype pS vs pN (by gene, by population)
ggsave(
  sub("\\.png$", ".pdf", OUT_WITHIN_MITO_PS_PN),
  p1,
  width = 18,
  height = 10,
  device = cairo_pdf
)

# 6) Within-population pS vs pN
ggsave(
  sub("\\.png$", ".pdf", OUT_WITHIN_POP_PS_PN),
  p2,
  width = 18,
  height = 10,
  device = cairo_pdf
)

# 7) Between-mitotype pairs pS vs pN (within population)
ggsave(
  sub("\\.png$", ".pdf", OUT_BETWEEN_PAIR_PS_PN),
  p3,
  width = 18,
  height = 10,
  device = cairo_pdf
)

# 8) Combined pS vs pN (all comparisons)
ggsave(
  sub("\\.png$", ".pdf", OUT_ALL_PS_PN_PNG),
  p_all_pspn,
  width = 16,
  height = 10,
  device = cairo_pdf
)

# 9) Within-mitotype pS vs pN (gene facets)
ggsave(
  sub("\\.png$", ".pdf", OUT_WITHIN_MITO_PS_PN_PNG),
  p_within_mito_pspn,
  width = 16,
  height = 10,
  device = cairo_pdf
)

# 10) Within-population pS vs pN (gene facets)
ggsave(
  sub("\\.png$", ".pdf", OUT_WITHIN_POP_PS_PN_PNG),
  p_within_pop_pspn,
  width = 16,
  height = 10,
  device = cairo_pdf
)

# 11) Between-mitotype pairs pS vs pN (gene facets)
ggsave(
  sub("\\.png$", ".pdf", OUT_BETWEEN_PAIRS_PS_PN_PNG),
  p_between_pairs_pspn,
  width = 16,
  height = 10,
  device = cairo_pdf
)

# 12) Pairwise all-mitotypes pS vs pN
ggsave(
  sub("\\.png$", ".pdf", OUT_PAIRWISE_ALL_PNG),
  p_pairwise_all,
  width = 16,
  height = 10,
  device = cairo_pdf
)