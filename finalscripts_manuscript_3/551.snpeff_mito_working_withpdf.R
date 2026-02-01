#ijob -A berglandlab -c10 -p standard --mem=40G
#module load gcc/11.4.0  openmpi/4.1.4 icu R/4.3.1;R


  library(vcfR)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(readr)

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

PSEUDO_LOC <- 0.1

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
# KEEP BIALLELIC SITES (ALT has no comma)
# ============================================================
bial_mask <- !is.na(fix$ALT) & !str_detect(fix$ALT, ",")
fix_bi <- fix[bial_mask, , drop = FALSE]
gt_bi  <- gt[bial_mask, , drop = FALSE]
log_msg("Biallelic records kept: ", nrow(fix_bi))

# ============================================================
# PARSE snpEff ANN -> gene + effect (syn/nonsyn)
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

    if (str_detect(eff, "missense_variant")) {
      return(c(gene, "nonsyn"))
    }
    if (str_detect(eff, "synonymous_variant") && is.na(best_eff)) {
      best_gene <- gene
      best_eff  <- "syn"
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
ggsave(sub("\\.png$", ".pdf", OUT_WITHIN_PNG), p_within, width = 14, height = 10, device = cairo_pdf)
log_msg("Wrote: ", OUT_WITHIN_PNG)

# ============================================================
# Within-location, within-mitotype PN/PS
# ============================================================
OUT_LOC_WITHIN_TSV <- file.path(OUT_DIR, paste0(PREFIX, ".within.location_mitotype.PNPS_by_gene.tsv"))
OUT_LOC_WITHIN_PNG <- file.path(OUT_DIR, paste0(PREFIX, ".within.location_mitotype.PNPS_by_gene.png"))

# ---- sample -> accuratelocation (standardize metadata clone key)

metadata_file  <- "/project/berglandlab/Robert/UKSequencing2022_2024/old_stuff/metadata_with_clone.csv"

metadata <- read.csv(metadata_file) %>%
  filter(!clone %in% c("Blank", "BLANK")) %>%
  mutate(
    Sample_ID = paste0(Well, "_clone"),
    accuratelocation = trimws(accuratelocation)
  ) %>%
  select(Sample_ID, accuratelocation)

meta_loc <- metadata %>%
  transmute(
    sample = str_remove(as.character(Sample_ID), "_clone$"),
    accuratelocation = na_if(as.character(accuratelocation), "")
  ) %>%
  distinct()

sample_map <- mito_df %>%
  mutate(sample = as.character(sample),
         mitotype = as.character(mitotype)) %>%
  left_join(meta_loc, by = "sample") %>%
  filter(!is.na(accuratelocation), accuratelocation != "")

keep_loc_samples <- intersect(colnames(allele_mat), sample_map$sample)
allele_loc <- allele_mat[, keep_loc_samples, drop = FALSE]
sample_map <- sample_map %>% filter(sample %in% keep_loc_samples)

log_msg("Samples with accuratelocation: ", nrow(sample_map),
        " across ", n_distinct(sample_map$accuratelocation), " locations.")

is_poly_within <- function(allele_mat, samples) {
  if (length(samples) < 2) return(rep(FALSE, nrow(allele_mat)))
  sub <- allele_mat[, samples, drop = FALSE]
  apply(sub, 1, function(v) {
    vv <- v[!is.na(v)]
    length(unique(vv)) > 1
  })
}

loc_levels   <- sort(unique(sample_map$accuratelocation))
mito_levels2 <- sort(unique(sample_map$mitotype))

hits <- vector("list", length(loc_levels) * length(mito_levels2))
ii <- 0L

for (loc in loc_levels) {
  loc_df <- sample_map %>% filter(accuratelocation == loc)

  for (m in mito_levels2) {
    s <- loc_df %>% filter(mitotype == m) %>% pull(sample)
    if (length(s) < 2) next

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
  arrange(accuratelocation, gene, mitotype) %>%
  mutate(label = paste0("pN=", nonsyn, ", pS=", syn))

write_tsv(within_loc_mito_pnps, OUT_LOC_WITHIN_TSV)
log_msg("Wrote: ", OUT_LOC_WITHIN_TSV)

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
ggsave(sub("\\.png$", ".pdf", OUT_LOC_WITHIN_PNG), p_loc_within, width = 18, height = 10, device = cairo_pdf)
log_msg("Wrote: ", OUT_LOC_WITHIN_PNG)

# ============================================================
# (B) WITHIN EACH POPULATION (location): PN/PS by gene
# ============================================================
OUT_POP_WITHIN_TSV <- file.path(OUT_DIR, paste0(PREFIX, ".within_population.PNPS_by_gene.tsv"))
OUT_POP_WITHIN_PNG <- file.path(OUT_DIR, paste0(PREFIX, ".within_population.PNPS_by_gene.png"))

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
ggsave(sub("\\.png$", ".pdf", OUT_POP_WITHIN_PNG), p_within_pop, width = 18, height = 10, device = cairo_pdf)
log_msg("Wrote: ", OUT_POP_WITHIN_PNG)

# ============================================================
# (C) BETWEEN MITOTYPES WITHIN EACH POPULATION: PN/PS by gene
# ============================================================
OUT_POP_BETWEEN_TSV <- file.path(OUT_DIR, paste0(PREFIX, ".between_mitotypes_within_population.PNPS_by_gene.tsv"))
OUT_POP_BETWEEN_PNG <- file.path(OUT_DIR, paste0(PREFIX, ".between_mitotypes_within_population.PNPS_by_gene.png"))

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
ggsave(sub("\\.png$", ".pdf", OUT_POP_BETWEEN_PNG), p_between_pop, width = 30, height = 10, device = cairo_pdf)
log_msg("Wrote: ", OUT_POP_BETWEEN_PNG)

# ============================================================
# Helper: plot pS vs pN (save PNG + PDF)
# ============================================================
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
  ggsave(sub("\\.png$", ".pdf", out_png), gg, width = 18, height = 10, device = cairo_pdf)
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

# ============================================================
# Combined pS vs pN (all comparisons) + PNG/PDF
# ============================================================
OUT_ALL_PS_PN_TSV <- file.path(OUT_DIR, paste0(PREFIX, ".all_comparisons.ps_vs_pn.tsv"))
OUT_ALL_PS_PN_PNG <- file.path(OUT_DIR, paste0(PREFIX, ".all_comparisons.ps_vs_pn.by_gene.png"))

df_pspn_all <- bind_rows(
  within_loc_mito_pnps %>%
    transmute(gene, accuratelocation, comparison = "Within mitotype", syn, nonsyn),
  within_pop_pnps %>%
    transmute(gene, accuratelocation, comparison = "Within population", syn, nonsyn),
  between_pop_pnps %>%
    transmute(gene, accuratelocation, comparison = "Between mitotypes (pairwise)", syn, nonsyn)
) %>%
  filter(!is.na(gene), !is.na(accuratelocation)) %>%
  mutate(
    comparison = factor(comparison,
      levels = c("Within population", "Within mitotype", "Between mitotypes (pairwise)")
    ),
    pond = factor(accuratelocation)
  ) %>%
  filter((syn + nonsyn) > 0)

write_tsv(df_pspn_all, OUT_ALL_PS_PN_TSV)
log_msg("Wrote: ", OUT_ALL_PS_PN_TSV)

p_all_pspn <- ggplot(df_pspn_all, aes(x = nonsyn, y = syn, color = pond, shape = comparison)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.6, color = "grey40") +
  geom_point(size = 2.6, alpha = 0.9, position = position_jitter(width = 0.05, height = 0.05)) +
  facet_wrap(~ gene, scales = "free") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 0), legend.position = "right") +
  labs(
    x = "pN (nonsyn count)",
    y = "pS (syn count)",
    color = "Pond",
    shape = "Comparison",
    title = "pS vs pN by gene",
    subtitle = "Dashed line shows pN = pS (slope = 1)"
  )

ggsave(OUT_ALL_PS_PN_PNG, p_all_pspn, width = 16, height = 10, dpi = 300)
ggsave(sub("\\.png$", ".pdf", OUT_ALL_PS_PN_PNG), p_all_pspn, width = 16, height = 10, device = cairo_pdf)
log_msg("Wrote: ", OUT_ALL_PS_PN_PNG)

# ============================================================
# Helper: 3 separate pS vs pN figures (save PNG + PDF)
# ============================================================
plot_pspn_by_gene <- function(df, out_png, title, subtitle = NULL) {
  df <- df %>%
    filter(!is.na(gene), !is.na(accuratelocation)) %>%
    mutate(pond = factor(accuratelocation)) %>%
    filter((syn + nonsyn) > 0)

  p <- ggplot(df, aes(x = nonsyn, y = syn, color = pond)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.6, color = "grey40") +
    geom_point(size = 2.6, alpha = 0.9, position = position_jitter(width = 0.05, height = 0.05)) +
    facet_wrap(~ gene, scales = "free") +
    theme_bw() +
    theme(panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 0), legend.position = "right") +
    labs(
      x = "pN (nonsyn count)",
      y = "pS (syn count)",
      color = "Pond",
      title = title,
      subtitle = subtitle
    )

  ggsave(out_png, p, width = 16, height = 10, dpi = 300)
  ggsave(sub("\\.png$", ".pdf", out_png), p, width = 16, height = 10, device = cairo_pdf)
  log_msg("Wrote: ", out_png)
  p
}

df_within_mito <- within_loc_mito_pnps %>% transmute(gene, accuratelocation, syn, nonsyn)

OUT_WITHIN_MITO_PS_PN_PNG <- file.path(OUT_DIR, paste0(PREFIX, ".within_mitotype.ps_vs_pn.by_gene.png"))

p_within_mito_pspn <- plot_pspn_by_gene(
  df_within_mito,
  OUT_WITHIN_MITO_PS_PN_PNG,
  title = "Within mitotype (per pond): pS vs pN by gene",
  subtitle = "Dashed line shows pN = pS (slope = 1)"
)

df_within_pop <- within_pop_pnps %>% transmute(gene, accuratelocation, syn, nonsyn)

OUT_WITHIN_POP_PS_PN_PNG <- file.path(OUT_DIR, paste0(PREFIX, ".within_population.ps_vs_pn.by_gene.png"))

p_within_pop_pspn <- plot_pspn_by_gene(
  df_within_pop,
  OUT_WITHIN_POP_PS_PN_PNG,
  title = "Within pond: pS vs pN by gene",
  subtitle = "Dashed line shows pN = pS (slope = 1)"
)

df_between_pairs <- between_pop_pnps %>%
  transmute(gene, accuratelocation, pair, syn, nonsyn) %>%
  filter((syn + nonsyn) > 0) %>%
  mutate(pond = factor(accuratelocation), pair = factor(pair))

OUT_BETWEEN_PAIRS_PS_PN_PNG <- file.path(OUT_DIR, paste0(PREFIX, ".between_mitotype_pairs.ps_vs_pn.by_gene.png"))

p_between_pairs_pspn <- ggplot(df_between_pairs, aes(x = nonsyn, y = syn, color = pond, shape = pair)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.6, color = "grey40") +
  geom_point(size = 2.6, alpha = 0.9, position = position_jitter(width = 0.05, height = 0.05)) +
  facet_wrap(~ gene, scales = "free") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), legend.position = "right") +
  xlim(0, 30) + ylim(0, 80) +
  labs(
    x = "pN (nonsyn count)",
    y = "pS (syn count)",
    color = "Pond",
    shape = "Mitotype pair",
    title = "Between mitotype pairs (within pond): pS vs pN by gene",
    subtitle = "Each point = one mitotype comparison within a pond; dashed line is pN = pS"
  )

ggsave(OUT_BETWEEN_PAIRS_PS_PN_PNG, p_between_pairs_pspn, width = 16, height = 10, dpi = 300)
ggsave(sub("\\.png$", ".pdf", OUT_BETWEEN_PAIRS_PS_PN_PNG), p_between_pairs_pspn, width = 16, height = 10, device = cairo_pdf)
log_msg("Wrote: ", OUT_BETWEEN_PAIRS_PS_PN_PNG)

# ============================================================
# Pairwise between all mitotypes (ignoring population): save PNG + PDF
# ============================================================
PSEUDO_PAIR <- if (exists("PSEUDO_LOC")) PSEUDO_LOC else 0.1

OUT_PAIRWISE_ALL_TSV <- file.path(OUT_DIR, paste0(PREFIX, ".pairwise_all_mitotypes.PNPS_by_gene.tsv"))
OUT_PAIRWISE_ALL_PNG <- file.path(OUT_DIR, paste0(PREFIX, ".pairwise_all_mitotypes.ps_vs_pn.by_gene.png"))

stopifnot(exists("mito_levels"), exists("fixed_by_mito"), exists("var_anno"))
stopifnot(all(c("variant_index","gene","effect") %in% names(var_anno)))
stopifnot(nrow(fixed_by_mito) == nrow(var_anno))

pair_mat <- t(combn(mito_levels, 2))
pair_df  <- tibble(mitoA = pair_mat[, 1], mitoB = pair_mat[, 2]) %>%
  mutate(pair = paste0(mitoA, " vs ", mitoB))

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

p_pairwise_all <- ggplot(pairwise_all, aes(x = pN, y = pS, color = pair)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.6, color = "grey40") +
  geom_point(size = 2.4, alpha = 0.85, position = position_jitter(width = 0.05, height = 0.05)) +
  facet_wrap(~ gene, scales = "free") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 0), legend.position = "right") +
  labs(
    x = "pN (nonsyn fixed-diff count)",
    y = "pS (syn fixed-diff count)",
    color = "Mitotype pair",
    title = "Pairwise between all mitotypes: pS vs pN by gene (ignoring population)",
    subtitle = "Each point = one mitotype pair; variants counted when both mitotypes are fixed and differ; dashed line is pN = pS"
  )

ggsave(OUT_PAIRWISE_ALL_PNG, p_pairwise_all, width = 16, height = 10, dpi = 300)
ggsave(sub("\\.png$", ".pdf", OUT_PAIRWISE_ALL_PNG), p_pairwise_all, width = 16, height = 10, device = cairo_pdf)
log_msg("Wrote: ", OUT_PAIRWISE_ALL_PNG)

# ============================================================
# Pairwise all-mitotypes, NA groups: save PNG + PDF
# ============================================================
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
    comparison_group = factor(comparison_group, levels = c("NA2", "NA1", "Between NA1 and NA2"))
  )

OUT_PAIRWISE_ALL_PNG <- file.path(
  OUT_DIR, paste0(PREFIX, ".pairwise_all_mitotypes.ps_vs_pn.by_gene.NA_groups.png")
)

p_pairwise_all <- ggplot(
  pairwise_all,
  aes(x = pN, y = pS, color = pair, shape = comparison_group)
) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.6, color = "grey40") +
  geom_point(size = 2.6, alpha = 0.85, position = position_jitter(width = 0.05, height = 0.05)) +
  facet_wrap(~ gene, scales = "free") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 0), legend.position = "right") +
  labs(
    x = "pN (nonsyn fixed-difference count)",
    y = "pS (syn fixed-difference count)",
    color = "Mitotype pair",
    shape = "Comparison group",
    title = "Pairwise pS vs pN between all mitotypes (by gene)",
    subtitle = "Shapes indicate NA groups; dashed line shows pN = pS"
  )

ggsave(OUT_PAIRWISE_ALL_PNG, p_pairwise_all, width = 16, height = 10, dpi = 300)
ggsave(sub("\\.png$", ".pdf", OUT_PAIRWISE_ALL_PNG), p_pairwise_all, width = 16, height = 10, device = cairo_pdf)
log_msg("Wrote: ", OUT_PAIRWISE_ALL_PNG)





































suppressPackageStartupMessages({
  library(vcfR)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(readr)
})

# -------------------------
# EDIT PATHS
# -------------------------
VCF_IN   <- "/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/snpeff/output.ann.vcf"
MITO_TSV <- "/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/mito_types.csv"

OUT_DIR  <- "/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/snpeff"
PREFIX   <- "biallelic.ABCF_vs_other"

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# -------------------------
# Settings
# -------------------------
GROUP1 <- c("A","B","C","F")  # merge these into one group
PSEUDO <- 0.5                 # pseudocount for PN/PS

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
allele_mat <- apply(gt2, c(1,2), gt_to_allele)

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
















suppressPackageStartupMessages({
  library(vcfR)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(readr)
})

# -------------------------
# EDIT PATHS
# -------------------------
VCF_IN   <- "/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/snpeff/output.ann.vcf"
MITO_TSV <- "/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/mito_types.csv"

OUT_DIR  <- "/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/snpeff"
PREFIX   <- "biallelic.pairwise.mergedGroups"

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# -------------------------
# Settings
# -------------------------
GROUP1 <- c("A","B","C","F")  # merged into ABCF
PSEUDO <- 0.5

OUT_PAIR_BETWEEN_TSV <- file.path(OUT_DIR, paste0(PREFIX, ".between_fixed.PNPS.tsv"))
OUT_PAIR_BETWEEN_PNG <- file.path(OUT_DIR, paste0(PREFIX, ".between_fixed.PNPS.png"))

OUT_PAIR_WITHIN_TSV  <- file.path(OUT_DIR, paste0(PREFIX, ".within_poly.PNPS.tsv"))
OUT_PAIR_WITHIN_PNG  <- file.path(OUT_DIR, paste0(PREFIX, ".within_poly.PNPS.png"))

# ============================================================
# (1) Load mitotype map: needs columns sample, mitotype
# ============================================================
mito_df_raw <- suppressMessages(readr::read_delim(
  MITO_TSV,
  delim = ifelse(grepl("\\.csv$", MITO_TSV, ignore.case = TRUE), ",", "\t"),
  show_col_types = FALSE
))

if (!all(c("sample","mitotype") %in% names(mito_df_raw))) {
  if (ncol(mito_df_raw) >= 2) {
    mito_df_raw <- mito_df_raw[, 1:2]
    colnames(mito_df_raw) <- c("sample","mitotype")
  } else stop("MITO_TSV must have at least 2 columns: sample, mitotype")
}

mito_df <- mito_df_raw %>%
  transmute(sample = as.character(sample),
            mitotype = as.character(mitotype)) %>%
  filter(!is.na(sample), !is.na(mitotype)) %>%
  distinct() %>%
  mutate(group = ifelse(mitotype %in% GROUP1, "ABCF", "Other"))

# ============================================================
# (2) Read VCF + keep BIALLELIC only
# ============================================================
vcf <- vcfR::read.vcfR(VCF_IN, verbose = FALSE)
gt_all <- vcfR::extract.gt(vcf, element = "GT")  # variants x samples
fix_all <- as.data.frame(vcf@fix, stringsAsFactors = FALSE)

bial_mask <- !is.na(fix_all$ALT) & !str_detect(fix_all$ALT, ",")
fix <- fix_all[bial_mask, , drop = FALSE]
gt  <- gt_all[bial_mask, , drop = FALSE]
message("Biallelic sites kept: ", nrow(fix))

# Keep only samples with mitotype labels
keep_samples <- intersect(colnames(gt), mito_df$sample)
if (length(keep_samples) < 2) stop("Too few VCF samples match mitotype map. Check sample naming.")
gt <- gt[, keep_samples, drop = FALSE]
mito_df <- mito_df %>% filter(sample %in% keep_samples)

mito_levels <- sort(unique(mito_df$mitotype))
mito_to_samples <- split(mito_df$sample, mito_df$mitotype)

# ============================================================
# (3) Parse ANN -> gene + syn/nonsyn
# ============================================================
parse_ann <- function(info) {
  if (is.na(info) || !str_detect(info, "ANN=")) return(c(NA_character_, NA_character_))
  ann <- str_replace(info, "^.*ANN=", "")
  ann <- str_split(ann, ";", simplify = TRUE)[1]
  entries <- str_split(ann, ",")[[1]]

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

info_vec <- vcfR::getINFO(vcf)[bial_mask]
anno_mat <- t(vapply(info_vec, parse_ann, FUN.VALUE = character(2)))
fix$gene   <- anno_mat[,1]
fix$effect <- anno_mat[,2]

fix2 <- fix %>% filter(effect %in% c("syn","nonsyn"), !is.na(gene))
if (nrow(fix2) == 0) stop("No syn/nonsyn annotations found in biallelic sites.")

row_pos_in_fix <- match(rownames(fix2), rownames(fix))
stopifnot(all(!is.na(row_pos_in_fix)))
gt2 <- gt[row_pos_in_fix, , drop = FALSE]
stopifnot(nrow(gt2) == nrow(fix2))

# Convert GT to allele index (take first allele if diploid encoding appears)
gt_to_allele <- function(x) {
  if (is.na(x) || x %in% c(".", "./.", ".|.")) return(NA_integer_)
  a <- str_split(x, "[/|]", simplify = TRUE)[1]
  suppressWarnings(as.integer(a))
}
allele_mat <- apply(gt2, c(1,2), gt_to_allele)

var_anno <- fix2 %>% transmute(variant_index = seq_len(nrow(fix2)), gene, effect)

# ============================================================
# (4A) Pairwise BETWEEN (fixed differences) between merged groups
# Here it's just ONE comparison: ABCF vs Other, but we still plot per gene.
# ============================================================
# fixed allele per merged group
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

p_between <- ggplot(pair_between, aes(x = pN, y = pS)) +
  geom_point(size = 2.6, alpha = 0.9) +
  facet_wrap(~ gene, scales = "free") +
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
        OUT_PAIR_BETWEEN_TSV, "\n", OUT_PAIR_BETWEEN_PNG)







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
    pN_pct = 100 * nonsyn / total,
    pS_pct = 100 * syn    / total,

    # optional PN/PS still retained
    pnps = (nonsyn + PSEUDO) / (syn + PSEUDO),
    log10_pnps = log10(pnps)
  ) %>%
  filter(total > 0) %>%
  arrange(gene)

readr::write_tsv(pair_between, OUT_PAIR_BETWEEN_TSV)

# -------------------------
# Plot: percent pN vs percent pS
# -------------------------
p_between <- ggplot(pair_between, aes(x = pN_pct, y = pS_pct, color = gene)) +
  geom_point(size = 2.6, alpha = 0.9) +
  coord_equal(xlim = c(0, 100), ylim = c(0, 100)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  labs(
    x = "pN (%) – nonsynonymous fixed differences",
    y = "pS (%) – synonymous fixed differences",
    color = "Gene",
    title = "Between merged groups (ABCF vs Other): pN% vs pS% by gene (biallelic)",
    subtitle = "Percentages within each gene; pN% + pS% = 100"
  )

ggsave(OUT_PAIR_BETWEEN_PNG, p_between, width = 14, height = 10, dpi = 300)
























suppressPackageStartupMessages({
  library(vcfR)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
})

# ---- paths (edit if needed) ----
VCF_IN   <- "/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/snpeff/output.ann.vcf"
MITO_TSV <- "/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/mito_types.csv"
OUT_DIR  <- "/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/snpeff"
PREFIX   <- "biallelic.snp_lists"

GROUP1 <- c("A","B","C","F")   # merged group ABCF

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

OUT_WITHIN_MITO_CSV   <- file.path(OUT_DIR, paste0(PREFIX, ".within_mitotypes.polymorphic.csv"))
OUT_PAIR_FIXED_CSV    <- file.path(OUT_DIR, paste0(PREFIX, ".between_mitotypes.fixed_differences.csv"))
OUT_GROUP_FIXED_CSV   <- file.path(OUT_DIR, paste0(PREFIX, ".between_groups.ABCF_vs_Other.fixed_differences.csv"))

# ============================================================
# (1) Load mitotype map (force first 2 cols if needed)
# ============================================================
mito_df_raw <- suppressMessages(read_delim(
  MITO_TSV,
  delim = ifelse(grepl("\\.csv$", MITO_TSV, ignore.case = TRUE), ",", "\t"),
  show_col_types = FALSE
))

if (!all(c("sample","mitotype") %in% names(mito_df_raw))) {
  if (ncol(mito_df_raw) < 2) stop("MITO_TSV must have >= 2 columns.")
  mito_df_raw <- mito_df_raw[, 1:2]
  names(mito_df_raw) <- c("sample","mitotype")
}

mito_df <- mito_df_raw %>%
  transmute(sample = as.character(sample),
            mitotype = as.character(mitotype)) %>%
  filter(!is.na(sample), !is.na(mitotype)) %>%
  distinct() %>%
  mutate(group = ifelse(mitotype %in% GROUP1, "ABCF", "Other"))

# ============================================================
# (2) Read VCF, subset to biallelic, subset samples
# ============================================================
vcf    <- read.vcfR(VCF_IN, verbose = FALSE)
gt_all <- extract.gt(vcf, element = "GT")                  # variants x samples
fix_all <- as.data.frame(vcf@fix, stringsAsFactors = FALSE)

bial_mask <- !is.na(fix_all$ALT) & !str_detect(fix_all$ALT, ",")
fix_bi <- fix_all[bial_mask, , drop = FALSE]
gt_bi  <- gt_all[bial_mask, , drop = FALSE]
info_bi <- getINFO(vcf)[bial_mask]

keep_samples <- intersect(colnames(gt_bi), mito_df$sample)
if (length(keep_samples) < 2) stop("Too few VCF samples match mitotype map.")
gt_bi <- gt_bi[, keep_samples, drop = FALSE]
mito_df <- mito_df %>% filter(sample %in% keep_samples)

# ============================================================
# (3) Parse ANN -> gene + effect (syn/nonsyn), prefer missense if present
# ============================================================
parse_ann <- function(info) {
  if (is.na(info) || !str_detect(info, "ANN=")) return(c(NA_character_, NA_character_))
  ann <- str_replace(info, "^.*ANN=", "")
  ann <- str_split(ann, ";", n = 2, simplify = TRUE)[1]
  entries <- str_split(ann, ",")[[1]]

  gene_syn <- NA_character_
  for (e in entries) {
    f <- str_split(e, "\\|", simplify = TRUE)
    if (ncol(f) < 5) next
    eff  <- f[1,2]
    gene <- f[1,4]
    if (is.na(gene) || gene == "") next

    if (str_detect(eff, "missense_variant")) return(c(gene, "nonsyn"))
    if (is.na(gene_syn) && str_detect(eff, "synonymous_variant")) gene_syn <- gene
  }
  if (!is.na(gene_syn)) return(c(gene_syn, "syn"))
  c(NA_character_, NA_character_)
}

anno_mat <- t(vapply(info_bi, parse_ann, FUN.VALUE = character(2)))
fix_bi$gene   <- anno_mat[,1]
fix_bi$effect <- anno_mat[,2]

keep_eff <- !is.na(fix_bi$gene) & fix_bi$effect %in% c("syn","nonsyn")
fix2 <- fix_bi[keep_eff, , drop = FALSE]
gt2  <- gt_bi[keep_eff, , drop = FALSE]

if (nrow(fix2) == 0) stop("No syn/nonsyn biallelic records after ANN parsing.")

# ============================================================
# (4) GT -> allele index (first allele)
# ============================================================
gt_to_allele <- function(x) {
  if (is.na(x) || x %in% c(".", "./.", ".|.")) return(NA_integer_)
  a <- str_split(x, "[/|]", simplify = TRUE)[1]
  suppressWarnings(as.integer(a))
}
allele_mat <- apply(gt2, c(1,2), gt_to_allele)

# ============================================================
# (5) Build a variant annotation table (one row per variant)
# ============================================================
var_tbl <- tibble(
  variant_index = seq_len(nrow(fix2)),
  CHROM = fix2$CHROM,
  POS   = as.integer(fix2$POS),
  REF   = fix2$REF,
  ALT   = fix2$ALT,
  variant_id = paste(fix2$CHROM, fix2$POS, fix2$REF, fix2$ALT, sep=":"),
  gene  = fix2$gene,
  effect = fix2$effect
)

# ============================================================
# (6) Helper: fixed allele and polymorphic flag within a sample set
# ============================================================
fixed_allele <- function(mat, samples) {
  sub <- mat[, samples, drop = FALSE]
  apply(sub, 1, function(v) {
    vv <- v[!is.na(v)]
    u <- unique(vv)
    if (length(u) == 1) u[1] else NA_integer_
  })
}

is_poly <- function(mat, samples) {
  if (length(samples) < 2) return(rep(FALSE, nrow(mat)))
  sub <- mat[, samples, drop = FALSE]
  apply(sub, 1, function(v) {
    vv <- v[!is.na(v)]
    length(unique(vv)) > 1
  })
}

# ============================================================
# (7) (A) SNPs polymorphic WITHIN each mitotype
# ============================================================
mito_levels <- sort(unique(mito_df$mitotype))
mito_to_samples <- split(mito_df$sample, mito_df$mitotype)

within_list <- lapply(mito_levels, function(m) {
  s <- mito_to_samples[[m]]
  if (length(s) < 2) return(NULL)
  idx <- which(is_poly(allele_mat, s))
  if (length(idx) == 0) return(NULL)

  var_tbl[idx, ] %>%
    mutate(mitotype = m) %>%
    select(mitotype, variant_id, CHROM, POS, REF, ALT, gene, effect)
})
within_mito_snps <- bind_rows(within_list)



write_csv(within_mito_snps, OUT_WITHIN_MITO_CSV)

# ============================================================
# (8) (B) SNPs FIXED DIFFERENCES between each mitotype pair
#     Condition: both mitotypes fixed (non-NA) and alleles differ
# ============================================================
fixed_by_mito <- sapply(mito_levels, function(m) fixed_allele(allele_mat, mito_to_samples[[m]]))
if (is.vector(fixed_by_mito)) fixed_by_mito <- matrix(fixed_by_mito, ncol = 1, dimnames = list(NULL, mito_levels))

pair_mat <- t(combn(mito_levels, 2))
pair_df <- tibble(mitoA = pair_mat[,1], mitoB = pair_mat[,2])

pair_rows <- vector("list", nrow(pair_df))
for (k in seq_len(nrow(pair_df))) {
  a <- pair_df$mitoA[k]; b <- pair_df$mitoB[k]
  fa <- fixed_by_mito[, a]
  fb <- fixed_by_mito[, b]
  idx <- which(!is.na(fa) & !is.na(fb) & fa != fb)
  if (length(idx) == 0) next

  pair_rows[[k]] <- var_tbl[idx, ] %>%
    mutate(
      mitoA = a, mitoB = b,
      pair = paste0(a, " vs ", b),
      alleleA = fa[idx],
      alleleB = fb[idx]
    ) %>%
    select(pair, mitoA, mitoB, variant_id, CHROM, POS, REF, ALT, gene, effect, alleleA, alleleB)
}

pair_fixed_snps <- bind_rows(pair_rows)
write_csv(pair_fixed_snps, OUT_PAIR_FIXED_CSV)

# ============================================================
# (9) (C) SNPs FIXED DIFFERENCES between merged groups: ABCF vs Other
# ============================================================
group_levels <- c("ABCF","Other")
group_to_samples <- split(mito_df$sample, mito_df$group)[group_levels]

fa_g <- fixed_allele(allele_mat, group_to_samples[["ABCF"]])
fb_g <- fixed_allele(allele_mat, group_to_samples[["Other"]])

idx_g <- which(!is.na(fa_g) & !is.na(fb_g) & fa_g != fb_g)

group_fixed_snps <- var_tbl[idx_g, ] %>%
  mutate(
    groupA = "ABCF",
    groupB = "Other",
    comparison = "ABCF vs Other",
    alleleA = fa_g[idx_g],
    alleleB = fb_g[idx_g]
  ) %>%
  select(comparison, groupA, groupB, variant_id, CHROM, POS, REF, ALT, gene, effect, alleleA, alleleB)

write_csv(group_fixed_snps, OUT_GROUP_FIXED_CSV)

message("Wrote CSVs:\n",
        OUT_WITHIN_MITO_CSV, "\n",
        OUT_PAIR_FIXED_CSV, "\n",
        OUT_GROUP_FIXED_CSV)
























OUT_MERGED_CSV <- file.path(OUT_DIR, paste0(PREFIX, ".ALL_SNP_LISTS.merged.csv"))
OUT_POS_PNG    <- file.path(OUT_DIR, paste0(PREFIX, ".ALL_SNP_LISTS.position.png"))
OUT_POS_PDF    <- sub("\\.png$", ".pdf", OUT_POS_PNG)

GROUP1 <- c("A","B","C","F") # ABCF
# If you already have mito_df$group from earlier, this will just re-derive consistently:
if (!("group" %in% names(mito_df))) {
  mito_df <- mito_df %>%
    mutate(group = if_else(mitotype %in% GROUP1, "ABCF", "Other"))
}

# ----------------------------
# Variant metadata table
# ----------------------------
var_meta <- fix2 %>%
  mutate(
    variant_index = seq_len(n()),
    POS = as.integer(POS),
    effect = as.character(effect),
    gene   = as.character(gene)
  ) %>%
  select(variant_index, CHROM, POS, REF, ALT, gene, effect)

stopifnot(nrow(var_meta) == nrow(allele_mat))

# ----------------------------
# Helper: fixed allele per set of samples
# ----------------------------
fixed_allele_vec <- function(allele_mat, samples) {
  if (length(samples) < 1) return(rep(NA_integer_, nrow(allele_mat)))
  sub <- allele_mat[, samples, drop = FALSE]
  apply(sub, 1, function(v) {
    vv <- v[!is.na(v)]
    u <- unique(vv)
    if (length(u) == 1) u[1] else NA_integer_
  })
}

# ----------------------------
# (1) Within each mitotype: polymorphic sites
# ----------------------------
mito_levels     <- sort(unique(mito_df$mitotype))
mito_to_samples <- split(mito_df$sample, mito_df$mitotype)

poly_by_mito <- matrix(FALSE, nrow = nrow(allele_mat), ncol = length(mito_levels),
                      dimnames = list(NULL, mito_levels))
for (m in mito_levels) {
  s <- mito_to_samples[[m]]
  if (length(s) < 2) next
  sub <- allele_mat[, s, drop = FALSE]
  poly_by_mito[, m] <- apply(sub, 1, function(v) {
    vv <- v[!is.na(v)]
    length(unique(vv)) > 1
  })
}

df_within_mito <- as.data.frame(poly_by_mito) %>%
  mutate(variant_index = seq_len(nrow(poly_by_mito))) %>%
  pivot_longer(-variant_index, names_to = "mitotype", values_to = "is_poly") %>%
  filter(is_poly) %>%
  transmute(
    variant_index,
    class = "within_mitotype_polymorphic",
    comparison = mitotype
  ) %>%
  left_join(var_meta, by = "variant_index")

# ----------------------------
# (2) Fixed SNPs between individual mitotypes (pairwise fixed differences)
#     Condition: both mitotypes fixed (no within-mito polymorphism among called),
#     and fixed alleles differ.
# ----------------------------
fixed_by_mito <- matrix(NA_integer_, nrow = nrow(allele_mat), ncol = length(mito_levels),
                        dimnames = list(NULL, mito_levels))
for (m in mito_levels) {
  s <- mito_to_samples[[m]]
  fixed_by_mito[, m] <- fixed_allele_vec(allele_mat, s)
}

pair_mat <- t(combn(mito_levels, 2))
pair_df  <- tibble(mitoA = pair_mat[,1], mitoB = pair_mat[,2]) %>%
  mutate(pair = paste0(mitoA, " vs ", mitoB))

pair_hits <- vector("list", nrow(pair_df))
for (k in seq_len(nrow(pair_df))) {
  a <- pair_df$mitoA[k]
  b <- pair_df$mitoB[k]
  fa <- fixed_by_mito[, a]
  fb <- fixed_by_mito[, b]
  idx <- which(!is.na(fa) & !is.na(fb) & fa != fb)
  if (length(idx) == 0) next

  pair_hits[[k]] <- tibble(
    variant_index = idx,
    class = "between_mitotypes_fixed",
    comparison = pair_df$pair[k]
  )
}

df_between_pairs <- bind_rows(pair_hits) %>%
  left_join(var_meta, by = "variant_index")

# ----------------------------
# (3) Fixed SNPs between ABCF vs Other (merged groups)
# ----------------------------
group_levels <- c("ABCF","Other")
group_to_samples <- split(mito_df$sample, mito_df$group)[group_levels]

fa_g <- fixed_allele_vec(allele_mat, group_to_samples[["ABCF"]])
fb_g <- fixed_allele_vec(allele_mat, group_to_samples[["Other"]])
idx_g <- which(!is.na(fa_g) & !is.na(fb_g) & fa_g != fb_g)

df_between_groups <- tibble(
  variant_index = idx_g,
  class = "between_groups_fixed",
  comparison = "ABCF vs Other"
) %>%
  left_join(var_meta, by = "variant_index")

# ----------------------------
# MERGE ALL INTO ONE CSV
# ----------------------------
all_snps <- bind_rows(df_within_mito, df_between_pairs, df_between_groups) %>%
  filter(!is.na(gene), effect %in% c("syn","nonsyn")) %>%
  mutate(
    effect = recode(effect, syn = "Synonymous", nonsyn = "Nonsynonymous"),
    class  = factor(class, levels = c("within_mitotype_polymorphic",
                                     "between_mitotypes_fixed",
                                     "between_groups_fixed"))
  ) %>%
  arrange(CHROM, POS, gene, class, comparison)

# write CSV
readr::write_csv(all_snps, OUT_MERGED_CSV)
message("Wrote merged CSV: ", OUT_MERGED_CSV)

# ============================================================
# PLOT: show where these differences sit along the genome
#   x = position
#   y = class (within vs between vs merged-group between)
#   color = effect (syn vs nonsyn)
#   faceted by gene (good for mito)
# ============================================================
p_pos <- ggplot(all_snps, aes(x = POS, y = class, color = effect)) +
  geom_point(alpha = 0.65, size = 1.4, position = position_jitter(height = 0.12, width = 0)) +
  facet_wrap(~ gene, scales = "free_x", ncol = 2) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 9)
  ) +
  labs(
    x = "Genomic position (POS)",
    y = "SNP class",
    color = "Effect",
    title = "Synonymous vs nonsynonymous SNPs by position",
    subtitle = "Merged table: within-mitotype polymorphism, pairwise fixed differences, and ABCF vs Other fixed differences (biallelic)"
  )

ggsave(OUT_POS_PNG, p_pos, width = 20, height = 15, dpi = 300)
ggsave(OUT_POS_PDF, p_pos, width = 20, height = 15, device = cairo_pdf)
message("Wrote position plot:\n", OUT_POS_PNG, "\n", OUT_POS_PDF)



between_groups <- (subset(all_snps, class == "between_groups_fixed"))


freq_one_var_dplyr <- between_groups %>%
  count(gene)
print(freq_one_var_dplyr)

# Frequency by two variables
freq_two_vars_dplyr <- between_groups %>%
  count(gene, effect)
print(freq_two_vars_dplyr)

freq_two_vars_dplyrnonsyn <- subset(freq_two_vars_dplyr, effect == "Nonsynonymous")

merged122 <- left_join(freq_one_var_dplyr, freq_two_vars_dplyrnonsyn, by = "gene")

merged122$n.y[is.na(merged122$n.y)] <- 0.1

merged122$nonsynovertotal <- merged122$n.y/merged122$n.x
merged122$synovertotal <- 1-(merged122$n.y/merged122$n.x)


p_pos_nonvssyn <- ggplot(merged122, aes(x = nonsynovertotal, y = synovertotal, color = gene)) +
  geom_point() +
  theme_bw() +
  labs(
    x = "prop nonsyn",
    y = "prop syn",
  ) + xlim(0,1)+ylim(0,1)

ggsave("/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/snpeff/nonsynoversyn.png", p_pos_nonvssyn, width = 9, height = 5, dpi = 300)