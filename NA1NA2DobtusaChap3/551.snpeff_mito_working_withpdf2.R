#ijob -A berglandlab -c10 -p standard --mem=40G
#module load gcc/11.4.0  openmpi/4.1.4 icu R/4.3.1;R







suppressPackageStartupMessages({
  library(vcfR)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(ggplot2)
})

# ---- paths (edit if needed) ----
VCF_IN   <- "/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/snpeff/output.ann.vcf"
MITO_TSV <- "/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/mito_types.csv"
OUT_DIR  <- "/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/snpeff"
PREFIX   <- "biallelic.snp_lists"

GROUP1 <- c("A","B","C")   # merged group ABC

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

OUT_WITHIN_MITO_CSV   <- file.path(OUT_DIR, paste0(PREFIX, ".within_mitotypes.polymorphic.csv"))
OUT_PAIR_FIXED_CSV    <- file.path(OUT_DIR, paste0(PREFIX, ".between_mitotypes.fixed_differences.csv"))
OUT_GROUP_FIXED_CSV   <- file.path(OUT_DIR, paste0(PREFIX, ".between_groups.ABC_vs_Other.fixed_differences.csv"))

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
  mutate(group = ifelse(mitotype %in% GROUP1, "ABC", "Other"))

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
# (9) (C) SNPs FIXED DIFFERENCES between merged groups: ABC vs Other
# ============================================================
group_levels <- c("ABC","Other")
group_to_samples <- split(mito_df$sample, mito_df$group)[group_levels]

fa_g <- fixed_allele(allele_mat, group_to_samples[["ABC"]])
fb_g <- fixed_allele(allele_mat, group_to_samples[["Other"]])

idx_g <- which(!is.na(fa_g) & !is.na(fb_g) & fa_g != fb_g)

group_fixed_snps <- var_tbl[idx_g, ] %>%
  mutate(
    groupA = "ABC",
    groupB = "Other",
    comparison = "ABC vs Other",
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

GROUP1 <- c("A","B","C") # ABC
# If you already have mito_df$group from earlier, this will just re-derive consistently:
if (!("group" %in% names(mito_df))) {
  mito_df <- mito_df %>%
    mutate(group = if_else(mitotype %in% GROUP1, "ABC", "Other"))
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
# (3) Fixed SNPs between ABC vs Other (merged groups)
# ----------------------------
group_levels <- c("ABC","Other")
group_to_samples <- split(mito_df$sample, mito_df$group)[group_levels]

fa_g <- fixed_allele_vec(allele_mat, group_to_samples[["ABC"]])
fb_g <- fixed_allele_vec(allele_mat, group_to_samples[["Other"]])
idx_g <- which(!is.na(fa_g) & !is.na(fb_g) & fa_g != fb_g)

df_between_groups <- tibble(
  variant_index = idx_g,
  class = "between_groups_fixed",
  comparison = "ABC vs Other"
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
    subtitle = "Merged table: within-mitotype polymorphism, pairwise fixed differences, and ABC vs Other fixed differences (biallelic)"
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

merged122$n.y[is.na(merged122$n.y)] <- 0

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


merged122$pNpS <- merged122$n.y / (merged122$n.x - merged122$n.y)
merged122$pN <- merged122$n.y 
merged122$pS <- (merged122$n.x - merged122$n.y)

merged122clean <- merged122

merged122clean$gene <- sub("^cds_transcript_", "", merged122clean$gene)


p_pos_nonvssyn2 <- ggplot(merged122clean, aes(x = gene, y = pNpS)) +
  geom_point() +
  geom_text(
    aes(label = paste0("N=", pN, "  S=", pS)),
    vjust = -0.6,
    size = 3,
    show.legend = FALSE
  ) +
  scale_x_discrete(expand = expansion(mult = c(0.05, 0.05))) +
  theme_bw() + ylim(0,1)+
  labs(
    x = "Gene",
    y = "pNpS"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  "/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/snpeff/ABCvOtherpNpS.png",
  p_pos_nonvssyn2,
  width = 11,
  height = 4,
  dpi = 300
)



ggsave(
  "/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/snpeff/ABCvOtherpNpS.pdf",
  p_pos_nonvssyn2,
  width = 11,
  height = 4,
  dpi = 300
)




















between_mitotypes_fixed <- (subset(all_snps, class == "between_mitotypes_fixed"))


freq_one_var_dplyr_mito <- between_mitotypes_fixed %>%
  count(gene, comparison)
print(freq_one_var_dplyr_mito)

# Frequency by two variables
freq_two_vars_dplyr_mito <- between_mitotypes_fixed %>%
  count(gene, effect, comparison)
print(freq_two_vars_dplyr_mito)

freq_two_vars_dplyr_mitononsyn <- subset(freq_two_vars_dplyr_mito, effect == "Nonsynonymous")
freq_two_vars_dplyr_mitosyn <- subset(freq_two_vars_dplyr_mito, effect == "Synonymous")

colnames(freq_two_vars_dplyr_mitononsyn) <- c("gene", "effect", "comparison", "numbernonsyn")
colnames(freq_two_vars_dplyr_mitosyn) <- c("gene", "effect", "comparison", "numbersyn")
colnames(freq_one_var_dplyr_mito) <- c("gene", "comparison", "total")


merged3 <- left_join(freq_one_var_dplyr_mito, freq_two_vars_dplyr_mitononsyn, by = c("gene", "comparison"))
merged4 <- left_join(merged3, freq_two_vars_dplyr_mitosyn, by = c("gene", "comparison"))

merged4$numbernonsyn[is.na(merged4$numbernonsyn)] <- 0
merged4$numbersyn[is.na(merged4$numbersyn)] <- 0.5

merged4$pNpS <- merged4$numbernonsyn/ merged4$numbersyn


merged4$gene <- sub("^cds_transcript_", "", merged4$gene)


p_pos_nonvssyn3 <- ggplot(merged4, aes(x = gene, y = pNpS)) +
  geom_point() + facet_wrap(~comparison)+
  geom_text(
    aes(label = paste0("N=", numbernonsyn, "  S=", numbersyn)),
    vjust = -0.6,
    size = 3,
    show.legend = FALSE
  ) +
  scale_x_discrete(expand = expansion(mult = c(0.05, 0.05))) +
  theme_bw() + ylim(0,5)+
  labs(
    x = "Gene",
    y = "pNpS"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  "/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/snpeff/ABCvOtherpNpS_comps.png",
  p_pos_nonvssyn3,
  width = 11,
  height = 11,
  dpi = 300
)

ggsave(
  "/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/snpeff/ABCvOtherpNpS.pdf",
  p_pos_nonvssyn3,
  width = 11,
  height = 4,
  dpi = 300
)








between_mitotypes_fixed <- (subset(all_snps, class == "between_mitotypes_fixed"))


freq_one_var_dplyr_mito <- between_mitotypes_fixed %>%
  count(comparison)
print(freq_one_var_dplyr_mito)

# Frequency by two variables
freq_two_vars_dplyr_mito <- between_mitotypes_fixed %>%
  count(effect, comparison)
print(freq_two_vars_dplyr_mito)

freq_two_vars_dplyr_mitononsyn <- subset(freq_two_vars_dplyr_mito, effect == "Nonsynonymous")
freq_two_vars_dplyr_mitosyn <- subset(freq_two_vars_dplyr_mito, effect == "Synonymous")

colnames(freq_two_vars_dplyr_mitononsyn) <- c("effect", "comparison", "numbernonsyn")
colnames(freq_two_vars_dplyr_mitosyn) <- c("effect", "comparison", "numbersyn")
colnames(freq_one_var_dplyr_mito) <- c("comparison", "total")


merged3 <- left_join(freq_one_var_dplyr_mito, freq_two_vars_dplyr_mitononsyn, by = c("comparison"))
merged4 <- left_join(merged3, freq_two_vars_dplyr_mitosyn, by = c("comparison"))

merged4$numbernonsyn[is.na(merged4$numbernonsyn)] <- 0
merged4$numbersyn[is.na(merged4$numbersyn)] <- 0.5

merged4$pNpS <- merged4$numbernonsyn/ merged4$numbersyn


merged4$gene <- sub("^cds_transcript_", "", merged4$gene)


p_pos_nonvssyn4 <- ggplot(merged4, aes(x = comparison, y = pNpS)) +
  geom_point()+
  geom_text(
    aes(label = paste0("N=", numbernonsyn, "  S=", numbersyn)),
    vjust = -0.6,
    size = 3,
    show.legend = FALSE
  ) +
  scale_x_discrete(expand = expansion(mult = c(0.05, 0.05))) +
  theme_bw() + ylim(0,0.5)+
  labs(
    x = "Gene",
    y = "pNpS"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  "/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/snpeff/ABCvOtherpNpS_comps_nogene.png",
  p_pos_nonvssyn4,
  width = 11,
  height = 11,
  dpi = 300
)



ggsave(
  "/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/snpeff/ABCvOtherpNpS_comps_nogene.pdf",
  p_pos_nonvssyn4,
  width = 11,
  height = 4,
  dpi = 300
)







p_pos_nonvssyn4 <- ggplot(merged4, aes(x = comparison, y = pNpS)) +
  geom_point()+
  geom_text(
    aes(label = paste0("N=", numbernonsyn, "  S=", numbersyn)),
    vjust = -0.6,
    size = 3,
    show.legend = FALSE
  ) +
  scale_x_discrete(expand = expansion(mult = c(0.05, 0.05))) +
  theme_bw() + ylim(0,5)+
  labs(
    x = "Gene",
    y = "pNpS"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  "/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/snpeff/ABCvOtherpNpS_comps_nogene.png",
  p_pos_nonvssyn4,
  width = 11,
  height = 11,
  dpi = 300
)










# classify comparisons into the 3 categories
merged4_class <- merged4 %>%
  mutate(
    a = str_extract(comparison, "^[A-Z]+"),
    b = str_extract(comparison, "[A-Z]+$"),
    group3 = case_when(
      a %in% c("A","B","C") & b %in% c("A","B","C") ~ "Within NA2 (A,B,C)",
      a %in% c("D","E","F") & b %in% c("D","E","F") ~ "Within NA1 (D,E,F)",
      (a %in% c("A","B","C") & b %in% c("D","E","F")) |
        (a %in% c("D","E","F") & b %in% c("A","B","C")) ~ "Between NA2 (ABC) and NA1 (DEF)",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(group3)) %>%
  mutate(group3 = factor(
    group3,
    levels = c("Within NA2 (A,B,C)", "Within NA1 (D,E,F)", "Between NA2 (ABC) and NA1 (DEF)")
  ))


library(patchwork)

 p_within_NA2 <-  ggplot(subset(merged4_class, group3 == "Within NA2 (A,B,C)"), aes(x = comparison, y = pNpS)) +
    geom_point() +
    geom_text(
      aes(label = paste0("N=", numbernonsyn, "  S=", numbersyn)),
      vjust = -0.6,
      size = 3,
      check_overlap = TRUE,
      show.legend = FALSE
    ) +
    scale_x_discrete(expand = expansion(mult = c(0.15, 0.15))) +
    coord_cartesian(ylim = c(0, 0.5)) +
    theme_bw() +
    labs(x = "Comparison", y = "pNpS") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(face = "bold")
    )  + ggtitle("Within NA2 (A,B,C)")


 p_within_NA1 <-  ggplot(subset(merged4_class, group3 == "Within NA1 (D,E,F)"), aes(x = comparison, y = pNpS)) +
    geom_point() +
    geom_text(
      aes(label = paste0("N=", numbernonsyn, "  S=", numbersyn)),
      vjust = -0.6,
      size = 3,
      check_overlap = TRUE,
      show.legend = FALSE
    ) +
    scale_x_discrete(expand = expansion(mult = c(0.15, 0.15))) +
    coord_cartesian(ylim = c(0, 0.5)) +
    theme_bw() +
    labs(x = "Comparison", y = "pNpS") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(face = "bold")
    )  + ggtitle("Within NA1 (D,E,F)")


 p_between <-  ggplot(subset(merged4_class, group3 == "Between NA2 (ABC) and NA1 (DEF)"), aes(x = comparison, y = pNpS)) +
    geom_point() +
    geom_text(
      aes(label = paste0("N=", numbernonsyn, "  S=", numbersyn)),
      vjust = -0.6,
      size = 3,
      check_overlap = TRUE,
      show.legend = FALSE
    ) +
    scale_x_discrete(expand = expansion(mult = c(0.08, 0.08))) +
    coord_cartesian(ylim = c(0, 0.5)) +
    theme_bw() +
    labs(x = "Comparison", y = "pNpS") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(face = "bold")
    ) + ggtitle("Between NA2 (ABC) and NA1 (DEF)")



p_pos_nonvssyn_faceted <- 
  (p_within_NA2 | p_within_NA1) /
  p_between +
  plot_layout(heights = c(1, 1.2)) +
  plot_annotation(tag_levels = "A")




ggsave(
  file.path("/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/snpeff/pNpS_comps_faceted.png"),
  p_pos_nonvssyn_faceted,
  width = 11,
  height = 6,
  dpi = 300
)

ggsave(
  file.path("/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/snpeff/pNpS_comps_faceted.pdf"),
  p_pos_nonvssyn_faceted,
  width = 11,
  height = 6
)
