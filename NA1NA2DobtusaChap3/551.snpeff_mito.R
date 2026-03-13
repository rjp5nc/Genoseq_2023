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

# -------------------------
# EDIT PATHS
# -------------------------
VCF_IN   <- "/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/snpeff/output.ann.vcf"
MITO_TSV <- "/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/mito_types.csv"  # or your mitotype map
OUT_TSV  <- "/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/snpeff/gene_syn_nonsyn_within_between.tsv"
OUT_PNG  <- "/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/snpeff/gene_syn_nonsyn_within_between.png"

# -------------------------
# Load mitotype map
# Expect columns: sample, mitotype
# If your file is clone_accuratelocation.csv etc, change mapping below.
# -------------------------
mito_df_raw <- suppressMessages(readr::read_delim(MITO_TSV, delim = ifelse(grepl("\\.csv$", MITO_TSV), ",", "\t"),
                                                 show_col_types = FALSE))
 


colnames(mito_df_raw) <- c("sample", "mitotype")

# ---- ADAPT THIS BLOCK to your file's column names ----
# If you already have a "sample" and "mitotype" column, you're done.
# Otherwise, set them here.
if (all(c("sample","mitotype") %in% names(mito_df_raw))) {
  mito_df <- mito_df_raw %>% select(sample, mitotype)
} else if (all(c("Well","clone") %in% names(mito_df_raw))) {
  # example: your earlier table had Well + clone; if "clone" is your sample ID and "Well" is mitotype, swap as needed.
  # CHANGE THESE TWO LINES to match your reality:
  mito_df <- mito_df_raw %>% transmute(sample = clone, mitotype = Well)
} else if (all(c("sample","haplotype") %in% names(mito_df_raw))) {
  mito_df <- mito_df_raw %>% transmute(sample = sample, mitotype = haplotype)
} else {
  stop("Couldn't find mitotype mapping columns. Make a 2-col file: sample<TAB>mitotype")
}

mito_df <- mito_df %>%
  filter(!is.na(sample), !is.na(mitotype)) %>%
  mutate(sample = as.character(sample),
         mitotype = as.character(mitotype)) %>%
  distinct()

# -------------------------
# Read VCF
# -------------------------
vcf <- vcfR::read.vcfR(VCF_IN, verbose = FALSE)
gt  <- vcfR::extract.gt(vcf, element = "GT")  # matrix: variants x samples
vcf_fix <- as.data.frame(vcf@fix, stringsAsFactors = FALSE) %>%
  mutate(variant_id = paste(CHROM, POS, REF, ALT, sep=":"))

samples <- colnames(gt)

# Keep only samples with mitotype labels
keep_samples <- intersect(samples, mito_df$sample)
if (length(keep_samples) < 2) stop("Too few VCF samples match mitotype map. Check sample naming.")

gt <- gt[, keep_samples, drop = FALSE]
mito_df <- mito_df %>% filter(sample %in% keep_samples)

# Helper: parse ANN -> choose gene + effect type
# ANN format: Allele|Annotation|Impact|Gene_Name|Gene_ID|Feature_Type|Feature_ID|...
# We'll:
#  - keep only entries containing missense_variant or synonymous_variant
#  - if multiple genes hit, pick the first Gene_Name (common for mito)
parse_ann_gene_effect <- function(info_vec) {
  # returns tibble(variant_row, gene, effect_class)
  # effect_class: "syn" or "nonsyn"
  out_gene <- rep(NA_character_, length(info_vec))
  out_eff  <- rep(NA_character_, length(info_vec))

  for (i in seq_along(info_vec)) {
    info <- info_vec[i]
    if (is.na(info) || !str_detect(info, "ANN=")) next

    ann_str <- str_replace(info, "^.*ANN=", "")
    ann_str <- str_split(ann_str, ";", n = 2, simplify = TRUE)[1] # up to first ;
    entries <- str_split(ann_str, ",")[[1]]

    # collect candidate rows
    cand <- lapply(entries, function(e) {
      fields <- str_split(e, "\\|", simplify = TRUE)
      if (ncol(fields) < 5) return(NULL)
      annot <- fields[1,2]
      gene  <- fields[1,4]
      if (is.na(annot) || is.na(gene) || gene == "") return(NULL)

      # annotation can be "missense_variant&something"
      parts <- unlist(str_split(annot, "&"))
      if ("missense_variant" %in% parts) return(list(gene=gene, eff="nonsyn"))
      if ("synonymous_variant" %in% parts) return(list(gene=gene, eff="syn"))
      NULL
    })
    cand <- Filter(Negate(is.null), cand)
    if (length(cand) == 0) next

    # prefer nonsyn if any
    effs <- vapply(cand, `[[`, "", "eff")
    genes <- vapply(cand, `[[`, "", "gene")
    if ("nonsyn" %in% effs) {
      out_eff[i] <- "nonsyn"
      out_gene[i] <- genes[which(effs=="nonsyn")[1]]
    } else {
      out_eff[i] <- "syn"
      out_gene[i] <- genes[1]
    }
  }

  tibble(gene = out_gene, effect = out_eff)
}

ann_tbl <- parse_ann_gene_effect(vcf_fix$INFO)
vcf_fix <- bind_cols(vcf_fix, ann_tbl)

# Keep only syn/nonsyn variants with gene name
vcf_fix2 <- vcf_fix %>%
  filter(!is.na(effect), !is.na(gene)) %>%
  mutate(row_id = row_number())

gt2 <- gt[vcf_fix2$row_id, , drop = FALSE]
stopifnot(nrow(gt2) == nrow(vcf_fix2))

# Convert GT strings to allele index (haploid):
# "0" ref, "1" alt, "./." missing, "." missing, "0/0" etc (just take first allele)
gt_to_allele <- function(x) {
  if (is.na(x) || x %in% c(".", "./.", ".|.")) return(NA_integer_)
  a <- str_split(x, "[/|]", simplify = TRUE)[1]
  suppressWarnings(as.integer(a))
}

allele_mat <- apply(gt2, c(1,2), gt_to_allele)

# Group sample indices by mitotype
mito_levels <- sort(unique(mito_df$mitotype))
mito_to_samples <- split(mito_df$sample, mito_df$mitotype)

# Determine within/between classification per variant
# within(m): polymorphic among called alleles within mitotype m
# between: each mitotype has a single fixed allele (among called), and at least two mitotypes have different fixed allele
within_long <- list()
between_flags <- logical(nrow(allele_mat))

fixed_by_mito <- matrix(NA_integer_, nrow=nrow(allele_mat), ncol=length(mito_levels),
                        dimnames=list(NULL, mito_levels))
poly_by_mito  <- matrix(FALSE, nrow=nrow(allele_mat), ncol=length(mito_levels),
                        dimnames=list(NULL, mito_levels))

for (j in seq_along(mito_levels)) {
  m <- mito_levels[j]
  s <- mito_to_samples[[m]]
  sub <- allele_mat[, s, drop = FALSE]
  # per row: set of called alleles
  aset <- apply(sub, 1, function(v) sort(unique(v[!is.na(v)])))
  poly <- vapply(aset, function(a) length(a) > 1, logical(1))
  fixed <- vapply(aset, function(a) if (length(a)==1) a[1] else NA_integer_, integer(1))
  poly_by_mito[, j]  <- poly
  fixed_by_mito[, j] <- fixed
}

# between: fixed within each mitotype (at least 2 mitotypes with non-NA fixed),
# and not all fixed alleles identical
between_flags <- apply(fixed_by_mito, 1, function(v) {
  vv <- v[!is.na(v)]
  if (length(vv) < 2) return(FALSE)
  length(unique(vv)) > 1
})

# Build tidy table of counts
# Within counts are per mitotype (a variant can contribute to within for multiple mitotypes)
within_df <- as.data.frame(poly_by_mito) %>%
  mutate(variant_index = seq_len(nrow(poly_by_mito))) %>%
  pivot_longer(-variant_index, names_to="mitotype", values_to="is_within") %>%
  filter(is_within)

between_df <- tibble(variant_index = which(between_flags)) %>%
  crossing(mitotype = mito_levels)  # to facet similarly; not strictly needed

# Join variant annotations
var_anno <- vcf_fix2 %>%
  transmute(variant_index = row_number(), gene, effect)

within_df2 <- within_df %>%
  left_join(var_anno, by="variant_index") %>%
  mutate(category = "within") %>%
  select(category, mitotype, gene, effect)

between_df2 <- between_df %>%
  left_join(var_anno, by="variant_index") %>%
  mutate(category = "between") %>%
  select(category, mitotype, gene, effect)

all_df <- bind_rows(within_df2, between_df2) %>%
  filter(!is.na(gene), !is.na(effect))

counts <- all_df %>%
  count(category, mitotype, gene, effect, name="n") %>%
  mutate(effect = recode(effect, syn="Synonymous", nonsyn="Nonsynonymous"),
         category = recode(category, within="Within mitotype", between="Between mitotypes"))

# write tidy counts
readr::write_tsv(counts, OUT_TSV)

# -------------------------
# Plot
# -------------------------
# Pick top genes by total to keep plot readable (mito genes are few, but just in case)
top_genes <- counts %>%
  group_by(gene) %>%
  summarise(total = sum(n), .groups="drop") %>%
  arrange(desc(total)) %>%
  slice_head(n = 20) %>%   # change if you want all genes
  pull(gene)

plot_df <- counts %>% filter(gene %in% top_genes)

p <- ggplot(plot_df, aes(x = gene, y = n, fill = effect)) +
  geom_col(position = "stack") +
  facet_grid(category ~ mitotype, scales = "free_y") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "Gene",
    y = "Variant count",
    fill = "Effect",
    title = "Synonymous vs nonsynonymous variants by gene",
    subtitle = "Within-mitotype polymorphism vs between-mitotype fixed differences"
  )

ggsave(OUT_PNG, p, width = 14, height = 8, dpi = 300)

cat("Wrote:\n", OUT_TSV, "\n", OUT_PNG, "\n", sep="")








stopifnot(all(c("gene","effect") %in% names(var_anno)))
stopifnot(nrow(fixed_by_mito) == nrow(var_anno))

# all unordered pairs of mitotypes
pair_df <- as.data.frame(t(combn(mito_levels, 2)), stringsAsFactors = FALSE)
colnames(pair_df) <- c("mitoA", "mitoB")

pair_counts <- vector("list", nrow(pair_df))

for (k in seq_len(nrow(pair_df))) {
  a <- pair_df$mitoA[k]
  b <- pair_df$mitoB[k]

  fa <- fixed_by_mito[, a]
  fb <- fixed_by_mito[, b]

  # variant is a between-difference for this pair if both fixed and different
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

# count syn/nonsyn per (pair,gene)
pair_summarised <- pair_long %>%
  count(pair, mitoA, mitoB, gene, effect, name = "n") %>%
  tidyr::pivot_wider(names_from = effect, values_from = n, values_fill = 0) %>%
  mutate(
    total = syn + nonsyn,
    prop_nonsyn = ifelse(total > 0, nonsyn / total, NA_real_),
    prop_syn    = ifelse(total > 0, syn / total, NA_real_)
  ) %>%
  filter(total > 0)

# optional: drop genes with very few pairwise differences to reduce noise
# pair_summarised <- pair_summarised %>% filter(total >= 5)



# Plot: one facet per gene, one point per pair
p_pair <- ggplot(pair_summarised, aes(x = prop_nonsyn, y = prop_syn, col = pair)) +
  geom_point(aes(size = total), alpha = 0.8) +
  # optional labels (can be cluttered)
  # ggrepel::geom_text_repel(aes(label = pair), size = 2.5, max.overlaps = 50)
  facet_wrap(~ gene) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "Proportion nonsynonymous (nonsyn / (syn + nonsyn))",
    y = "Proportion synonymous (syn / (syn + nonsyn))",
    size = "N (fixed differences)",
    title = "Pairwise between-mitotype fixed differences by gene",
    subtitle = "Each point = one mitotype pair; proportions within that pair and gene"
  ) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1))

# Save (single big faceted figure)
ggsave(OUT_PNG, p_pair, width = 14, height = 10, dpi = 300)

# Also save the underlying pairwise table if you want
readr::write_tsv(pair_summarised, OUT_TSV)

cat("Wrote:\n", OUT_TSV, "\n", OUT_PNG, "\n", sep="")


within_df <- as.data.frame(poly_by_mito) %>%
  mutate(variant_index = seq_len(nrow(poly_by_mito))) %>%
  pivot_longer(-variant_index, names_to="mitotype", values_to="is_within") %>%
  filter(is_within) %>%
  select(variant_index, mitotype)

var_anno <- fix2 %>%
  transmute(variant_index = seq_len(nrow(fix2)),
            gene = gene,
            effect = effect)

within_anno <- within_df %>%
  left_join(var_anno, by="variant_index") %>%
  filter(!is.na(gene), !is.na(effect))

# -------------------------
# Count syn/nonsyn within each mitotype by gene
# -------------------------
counts <- within_anno %>%
  mutate(effect = recode(effect, syn="Synonymous", nonsyn="Nonsynonymous")) %>%
  count(mitotype, gene, effect, name="n") %>%
  complete(mitotype, gene, effect, fill = list(n=0)) %>%
  arrange(mitotype, gene, effect)

readr::write_tsv(counts, OUT_TSV)

# -------------------------
# Plot (stacked bar by gene; facet by mitotype)
# -------------------------
# optional: keep top genes overall so plot stays readable
top_genes <- counts %>%
  group_by(gene) %>%
  summarise(total = sum(n), .groups="drop") %>%
  arrange(desc(total)) %>%
  slice_head(n = 20) %>%
  pull(gene)

plot_df <- counts %>% filter(gene %in% top_genes)

p <- ggplot(plot_df, aes(x = gene, y = n, fill = effect)) +
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
    title = "Within-mitotype polymorphism: synonymous vs nonsynonymous (by gene)"
  )

ggsave(OUT_PNG, p, width = 14, height = 8, dpi = 300)

cat("Wrote:\n", OUT_TSV, "\n", OUT_PNG, "\n", sep="")





p <- ggplot(plot_df, aes(x = gene, y = n, fill = effect)) +
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
    title = "Within-mitotype polymorphism: synonymous vs nonsynonymous (by gene)"
  )

ggsave(OUT_PNG, p, width = 14, height = 8, dpi = 300)







counts <- within_anno %>%
  count(mitotype, gene, effect, name="n") %>%
  complete(mitotype, gene, effect, fill = list(n=0)) %>%
  tidyr::pivot_wider(names_from = effect, values_from = n, values_fill = 0) %>%
  # effect columns are "syn" and "nonsyn" if present
  mutate(
    syn    = ifelse(is.na(syn), 0L, syn),
    nonsyn = ifelse(is.na(nonsyn), 0L, nonsyn),
    total = syn + nonsyn,
    prop_nonsyn = ifelse(total > 0, nonsyn / total, NA_real_),
    prop_syn    = ifelse(total > 0, syn    / total, NA_real_)
  ) %>%
  arrange(gene, mitotype)

readr::write_tsv(counts, OUT_TSV)

# -------------------------
# Plot: within-mitotype proportions (one point per mitotype per gene)
# -------------------------
plot_df <- counts %>% filter(total > 0)

p <- ggplot(plot_df, aes(x = prop_nonsyn, y = prop_syn, color = mitotype)) +
  geom_point(size = 2, alpha = 0.9) +
  coord_equal(xlim = c(0,1), ylim = c(0,1)) +
  facet_wrap(~ gene) +
  theme_bw() + xlim(0,1) + ylim(0,1)+
  theme(panel.grid.minor = element_blank()) +
  labs(
    x = "Proportion nonsynonymous (within mitotype)",
    y = "Proportion synonymous (within mitotype)",
    color = "Mitotype",
    title = "Within-mitotype syn vs nonsyn proportions by gene",
  )

ggsave(OUT_PNG, p, width = 14, height = 10, dpi = 300)


