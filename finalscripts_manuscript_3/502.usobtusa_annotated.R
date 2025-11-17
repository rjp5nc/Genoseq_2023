#ijob -A berglandlab -c10 -p standard --mem=80G
#module load gcc/11.4.0  openmpi/4.1.4 icu R/4.3.1
#R

library(SeqArray)
library(dplyr)
library(ggplot2)
library(tidyr)

# ------------------------------------------------------------
# 0) Inputs
# ------------------------------------------------------------

seqClose(g)

gds_file    <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/trimmed10bp_allsites_usobtusa_renamed_annotated.gds"
contig_file <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/contigmap.txt"

include_regions <- data.frame(
  contig = c(
    "JAACYE010000002.1", 
    "JAACYE010000003.1", "JAACYE010000003.1", "JAACYE010000003.1",
 #   "JAACYE010000006.1", "JAACYE010000006.1",
    "JAACYE010000007.1", "JAACYE010000007.1",
    "JAACYE010000008.1", "JAACYE010000008.1",
    "JAACYE010000011.1", "JAACYE010000011.1", "JAACYE010000011.1",
    "JAACYE010000012.1"
  ),

  #need the bit on chr 6  Rockpool4_F7
  win_start = c(
  4300000,
  2500000, 2600000, 7800000,
  1100000, 1200000, 
  1800000, 1900000, 
  1400000, 1500000, 6500000,
  1100000),
  win_end   = c(
  4399999,
  2599999, 2699999, 7899999,
  1199999, 1299999, 
  1899999, 1999999, 
  1499999, 1599999, 6599999,
  1199999)
)

include_regions <- include_regions %>%
  mutate(
    contig = as.character(contig),
    win_start = as.integer(win_start),
    win_end   = as.integer(win_end)
  )

include_regions$contig_mapped <- ifelse(
  include_regions$contig %in% contigmap$old,
  map_vec[ include_regions$contig ],
  include_regions$contig
)

include_regions


# g <- seqOpen(gds_file)

# chr_all <- seqGetData(g, "chromosome")
# pos_all <- seqGetData(g, "position")
# n_var   <- length(pos_all)

# # ------------------------------------------------------------
# # 2) Load contig map and map include_regions contigs
# # ------------------------------------------------------------

# contigmap <- read.table(contig_file, header = FALSE, stringsAsFactors = FALSE)
# colnames(contigmap) <- c("old", "new")
# map_vec <- setNames(contigmap$new, contigmap$old)

# include_regions <- include_regions %>%
#   mutate(
#     chr = ifelse(
#       contig %in% contigmap$old,
#       map_vec[contig],
#       contig
#     )
#   )

# # ------------------------------------------------------------
# # 3) Mark variants that fall in include_regions
# # ------------------------------------------------------------

# keep_var  <- rep(FALSE, n_var)
# region_id <- rep(NA_integer_, n_var)

# for (i in seq_len(nrow(include_regions))) {
#   ct <- include_regions$chr[i]
#   s  <- include_regions$win_start[i]
#   e  <- include_regions$win_end[i]

#   idx <- which(chr_all == ct & pos_all >= s & pos_all <= e)
#   if (length(idx)) {
#     keep_var[idx]  <- TRUE
#     region_id[idx] <- i
#   }
# }

# cat("Variants in include_regions:", sum(keep_var), "\n")

# # ------------------------------------------------------------
# # 4) Rebuild per-variant ANN strings
# # ------------------------------------------------------------

# ann_raw <- seqGetData(g, "annotation/info/ANN")
# len <- ann_raw$length
# dat <- ann_raw$data

# stopifnot(length(len) == n_var)

# idx_end   <- cumsum(len)
# idx_start <- c(1, head(idx_end, -1) + 1)

# ann_vec <- character(n_var)
# for (i in seq_along(len)) {
#   if (len[i] == 0L) {
#     ann_vec[i] <- NA_character_
#   } else {
#     ann_vec[i] <- paste(dat[idx_start[i]:idx_end[i]], collapse = ",")
#   }
# }

# has_ann <- !is.na(ann_vec)
# cat("Variants with any ANN:", sum(has_ann), "\n")

# in_region_with_ann <- keep_var & has_ann
# cat("Annotated variants in include_regions:", sum(in_region_with_ann), "\n")

# if (sum(in_region_with_ann) == 0) {
#   stop("No annotated variants in the specified include_regions.")
# }

# idx_keep <- which(in_region_with_ann)

# chr_sub       <- chr_all[idx_keep]
# pos_sub       <- pos_all[idx_keep]
# ann_sub       <- ann_vec[idx_keep]
# region_id_sub <- region_id[idx_keep]

# # ------------------------------------------------------------
# # 5) Helpers to extract effects and biotypes from ANN
# # ------------------------------------------------------------

# extract_effects <- function(ann_string) {
#   if (is.na(ann_string) || ann_string == "") return(NA_character_)
#   entries <- strsplit(ann_string, ",", fixed = TRUE)[[1]]
#   effects <- vapply(entries, function(e) {
#     fields <- strsplit(e, "\\|")[[1]]
#     if (length(fields) >= 2) fields[2] else NA_character_
#   }, character(1))
#   unname(effects)
# }

# extract_biotypes <- function(ann_string) {
#   if (is.na(ann_string) || ann_string == "") return(NA_character_)
#   entries <- strsplit(ann_string, ",", fixed = TRUE)[[1]]
#   biotypes <- vapply(entries, function(e) {
#     fields <- strsplit(e, "\\|")[[1]]
#     if (length(fields) >= 8) fields[8] else NA_character_
#   }, character(1))
#   unname(biotypes)
# }

# # primary effect = first non-NA effect
# ann_type <- vapply(ann_sub, function(x) {
#   eff <- extract_effects(x)
#   if (all(is.na(eff))) NA_character_ else eff[!is.na(eff)][1]
# }, character(1))

# # primary biotype = first non-NA biotype
# biotype <- vapply(ann_sub, function(x) {
#   bt <- extract_biotypes(x)
#   if (all(is.na(bt))) NA_character_ else bt[!is.na(bt)][1]
# }, character(1))

# cat("Top annotation types:\n")
# print(head(sort(table(ann_type), decreasing = TRUE), 15))

# cat("Top biotypes:\n")
# print(head(sort(table(biotype), decreasing = TRUE), 15))

# # ------------------------------------------------------------
# # 6) Collapse ann_type into broader effect categories
# # ------------------------------------------------------------

# ann_cat <- case_when(
#   grepl("stop_gained",            ann_type) ~ "stop_gained",
#   grepl("stop_lost",              ann_type) ~ "stop_lost",
#   grepl("start_lost",             ann_type) ~ "start_lost",
#   grepl("frameshift_variant",     ann_type) ~ "frameshift",
#   grepl("missense_variant",       ann_type) ~ "missense",
#   grepl("synonymous_variant",     ann_type) ~ "synonymous",
#   grepl("inframe_deletion",       ann_type) ~ "inframe_del",
#   grepl("inframe_insertion",      ann_type) ~ "inframe_ins",
#   grepl("splice_acceptor_variant|splice_donor_variant", ann_type) ~ "splice_site",
#   grepl("splice_region_variant",  ann_type) ~ "splice_region",
#   grepl("exon_loss_variant",      ann_type) ~ "exon_loss",
#   grepl("downstream_gene_variant",ann_type) ~ "downstream",
#   grepl("upstream_gene_variant",  ann_type) ~ "upstream",
#   grepl("intergenic_region",      ann_type) ~ "intergenic",
#   grepl("intron_variant",         ann_type) ~ "intron",
#   TRUE                                         ~ "other"
# )

# ann_cat <- factor(
#   ann_cat,
#   levels = c(
#     "stop_gained", "stop_lost", "start_lost",
#     "frameshift", "missense", "synonymous",
#     "inframe_del", "inframe_ins",
#     "splice_site", "splice_region",
#     "exon_loss",
#     "downstream", "upstream",
#     "intron", "intergenic",
#     "other"
#   )
# )

# # ------------------------------------------------------------
# # 7) Build plotting data frame with region labels
# # ------------------------------------------------------------

# plot_df <- data.frame(
#   chr       = chr_sub,
#   pos       = pos_sub,
#   ann_type  = ann_type,
#   ann_cat   = ann_cat,
#   biotype   = biotype,
#   region_id = region_id_sub
# )

# region_labels <- include_regions %>%
#   mutate(region_id = seq_len(n())) %>%
#   mutate(label = paste0(chr, ":", win_start, "-", win_end)) %>%
#   select(region_id, label)

# plot_df <- left_join(plot_df, region_labels, by = "region_id")

# head(plot_df)

# # ------------------------------------------------------------
# # 8) Plot effect categories across include_regions
# # ------------------------------------------------------------

# p_effects <- ggplot(plot_df, aes(x = pos, y = 0, color = ann_cat)) +
#   geom_point(alpha = 0.8, position = position_jitter(height = 0.1)) +
#   facet_wrap(~ label, scales = "free_x", ncol = 1) +
#   theme_bw() +
#   theme(
#     axis.text.y  = element_blank(),
#     axis.ticks.y = element_blank(),
#     axis.title.y = element_blank()
#   ) +
#   labs(
#     title = "SnpEff effect types across selected regions",
#     x     = "Position (bp)",
#     color = "Effect category"
#   )

# ggsave("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/p_effects.pdf", plot = p_effects, width = 12, height = 8, dpi = 300)

# # ------------------------------------------------------------
# # 9) Plot gene biotypes across include_regions
# # ------------------------------------------------------------

# p_biotype <- ggplot(plot_df, aes(x = pos, y = 0, color = biotype)) +
#   geom_point(alpha = 0.8, position = position_jitter(height = 0.1)) +
#   facet_wrap(~ label, scales = "free_x", ncol = 1) +
#   theme_bw() +
#   theme(
#     axis.text.y  = element_blank(),
#     axis.ticks.y = element_blank(),
#     axis.title.y = element_blank()
#   ) +
#   labs(
#     title = "Gene biotypes across selected regions",
#     x     = "Position (bp)",
#     color = "Biotype"
#   )

# ggsave("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/p_biotype.pdf", plot = p_biotype, width = 12, height = 8, dpi = 300)



# extract_impacts <- function(ann_string) {
#   if (is.na(ann_string) || ann_string == "") return(NA_character_)
#   entries <- strsplit(ann_string, ",", fixed = TRUE)[[1]]
#   impacts <- vapply(entries, function(e) {
#     fields <- strsplit(e, "\\|")[[1]]
#     if (length(fields) >= 3) fields[3] else NA_character_
#   }, character(1))
#   unique(na.omit(impacts))
# }

# impact_primary <- vapply(ann_sub, function(x) {
#   im <- extract_impacts(x)
#   if (length(im) == 0) NA_character_ else im[1]
# }, character(1))

# plot_df$impact <- impact_primary



# impactplot <- ggplot(plot_df, aes(x = pos, y = 0, color = impact)) +
#   geom_point(position = position_jitter(height = 0.1), alpha = 0.8) +
#   facet_wrap(~ label, scales = "free_x", ncol = 1) +
#   theme_bw() +
#   theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
#         axis.title.y = element_blank()) +
#   labs(title = "Variant impact across selected regions",
#        x = "Position (bp)", color = "Impact")


# ggsave("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/impactplot.pdf", plot = impactplot, width = 12, height = 8, dpi = 300)




# effect_impact_counts <- plot_df %>%
#   count(label, ann_cat, impact)

# effect_impact_countsplot <- ggplot(effect_impact_counts, aes(x = label, y = n, fill = ann_cat)) +
#   geom_col(position = "stack") +
#   coord_flip() +
#   theme_bw() +
#   labs(title = "Effect categories per region",
#        x = "Region", y = "Variant count", fill = "Effect")

# ggsave("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/effect_impact_countsplot.pdf", plot = effect_impact_countsplot, width = 12, height = 8, dpi = 300)











# gt <- seqGetData(g, "genotype")  # [ploidy, sample, variant]


# het_prop <- apply(gt, 3, function(G) {

#   a1 <- G[1, ]
#   a2 <- G[2, ]

#   # TRUE where both alleles are non-missing (≥ 0)
#   nonmiss <- (!is.na(a1)) & (!is.na(a2)) & (a1 >= 0) & (a2 >= 0)

#   # if nothing left, return NA
#   if (!any(nonmiss, na.rm = TRUE))
#     return(NA_real_)

#   het <- (a1 != a2) & nonmiss
#   sum(het, na.rm = TRUE) / sum(nonmiss, na.rm = TRUE)
# })
























gds_file    <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/trimmed10bp_allsites_usobtusa_renamed_annotated.gds"
contig_file <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/contigmap.txt"
out_dir     <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv"

## include_regions should already exist, e.g.:
## include_regions <- data.frame(contig = ..., win_start = ..., win_end = ...)

genomic_types <- read.csv("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/genomic_types_v2.csv")
genomic_types <- subset(genomic_types, Group == "A")

signatures <- read.csv("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/Gilmer_A_het_blocks_signature.csv")


merged_signatures <- genomic_types %>%
  left_join(
    signatures %>%
      select(sample, signature),
    by = c("CloneA" = "sample")
  )

genomic_types <- subset(merged_signatures, signature == 1111101)

## ============================================================
## 1) Open GDS and basic variant info
## ============================================================

g <- seqOpen(gds_file)

chr_all <- seqGetData(g, "chromosome")
pos_all <- seqGetData(g, "position")
n_var   <- length(pos_all)

samples_all <- seqGetData(g, "sample.id")
n_samp_all  <- length(samples_all)

cat("Total variants in GDS:", n_var, "\n")
cat("Total samples in GDS :", n_samp_all, "\n")

## ============================================================
## 2) Load contig map and map include_regions contigs
## ============================================================

contigmap <- read.table(contig_file, header = FALSE, stringsAsFactors = FALSE)
colnames(contigmap) <- c("old", "new")
map_vec <- setNames(contigmap$new, contigmap$old)

include_regions <- include_regions %>%
  mutate(
    chr = ifelse(
      contig %in% contigmap$old,
      map_vec[contig],
      contig
    )
  )

cat("include_regions:\n")
print(include_regions)

## ============================================================
## 3) Mark variants that fall in include_regions
## ============================================================

keep_var  <- rep(FALSE, n_var)
region_id <- rep(NA_integer_, n_var)

for (i in seq_len(nrow(include_regions))) {
  ct <- include_regions$chr[i]
  s  <- include_regions$win_start[i]
  e  <- include_regions$win_end[i]

  idx <- which(chr_all == ct & pos_all >= s & pos_all <= e)
  if (length(idx)) {
    keep_var[idx]  <- TRUE
    region_id[idx] <- i
  }
}

cat("Variants in include_regions:", sum(keep_var), "\n")

## ============================================================
## 4) Rebuild SnpEff ANN strings per variant
## ============================================================

ann_raw <- seqGetData(g, "annotation/info/ANN")
len <- ann_raw$length
dat <- ann_raw$data

stopifnot(length(len) == n_var)

idx_end   <- cumsum(len)
idx_start <- c(1, head(idx_end, -1) + 1)

ann_vec <- character(n_var)
for (i in seq_along(len)) {
  if (len[i] == 0L) {
    ann_vec[i] <- NA_character_
  } else {
    ann_vec[i] <- paste(dat[idx_start[i]:idx_end[i]], collapse = ",")
  }
}

has_ann <- !is.na(ann_vec)
cat("Variants with any ANN:", sum(has_ann), "\n")

in_region_with_ann <- keep_var & has_ann
cat("Annotated variants in include_regions:", sum(in_region_with_ann), "\n")

if (sum(in_region_with_ann) == 0) {
  stop("No annotated variants in the specified include_regions.")
}

idx_keep <- which(in_region_with_ann)

chr_sub       <- chr_all[idx_keep]
pos_sub       <- pos_all[idx_keep]
ann_sub       <- ann_vec[idx_keep]
region_id_sub <- region_id[idx_keep]

## ============================================================
## 5) Helper functions to parse ANN
## ============================================================

extract_effects <- function(ann_string) {
  if (is.na(ann_string) || ann_string == "") return(NA_character_)
  entries <- strsplit(ann_string, ",", fixed = TRUE)[[1]]
  effects <- vapply(entries, function(e) {
    fields <- strsplit(e, "\\|")[[1]]
    if (length(fields) >= 2) fields[2] else NA_character_
  }, character(1))
  unname(effects)
}

extract_impacts <- function(ann_string) {
  if (is.na(ann_string) || ann_string == "") return(NA_character_)
  entries <- strsplit(ann_string, ",", fixed = TRUE)[[1]]
  impacts <- vapply(entries, function(e) {
    fields <- strsplit(e, "\\|")[[1]]
    if (length(fields) >= 3) fields[3] else NA_character_
  }, character(1))
  unique(na.omit(unname(impacts)))
}

extract_biotypes <- function(ann_string) {
  if (is.na(ann_string) || ann_string == "") return(NA_character_)
  entries <- strsplit(ann_string, ",", fixed = TRUE)[[1]]
  biotypes <- vapply(entries, function(e) {
    fields <- strsplit(e, "\\|")[[1]]
    if (length(fields) >= 8) fields[8] else NA_character_
  }, character(1))
  unique(na.omit(unname(biotypes)))
}

extract_gene_ids <- function(ann_string) {
  if (is.na(ann_string) || ann_string == "") return(NA_character_)
  entries <- strsplit(ann_string, ",", fixed = TRUE)[[1]]
  gene_ids <- vapply(entries, function(e) {
    fields <- strsplit(e, "\\|")[[1]]
    if (length(fields) >= 5) fields[5] else NA_character_
  }, character(1))
  unique(na.omit(unname(gene_ids)))
}

## ============================================================
## 6) Primary effect, impact, biotype, genes per variant
## ============================================================

ann_type <- vapply(ann_sub, function(x) {
  eff <- extract_effects(x)
  if (all(is.na(eff))) NA_character_ else eff[!is.na(eff)][1]
}, character(1))

impact <- vapply(ann_sub, function(x) {
  im <- extract_impacts(x)
  if (length(im) == 0) NA_character_ else im[1]
}, character(1))

biotype <- vapply(ann_sub, function(x) {
  bt <- extract_biotypes(x)
  if (length(bt) == 0) NA_character_ else bt[1]
}, character(1))

gene_list_by_variant <- lapply(ann_sub, extract_gene_ids)
genes_in_regions     <- unique(unlist(gene_list_by_variant))

cat("Number of unique genes in include_regions (any samples):", length(genes_in_regions), "\n")
cat("First few genes:\n")
print(head(genes_in_regions))

## ============================================================
## 7) Collapse effects into broader categories
## ============================================================

ann_cat <- case_when(
  grepl("stop_gained",            ann_type) ~ "stop_gained",
  grepl("stop_lost",              ann_type) ~ "stop_lost",
  grepl("start_lost",             ann_type) ~ "start_lost",
  grepl("frameshift_variant",     ann_type) ~ "frameshift",
  grepl("missense_variant",       ann_type) ~ "missense",
  grepl("synonymous_variant",     ann_type) ~ "synonymous",
  grepl("inframe_deletion",       ann_type) ~ "inframe_del",
  grepl("inframe_insertion",      ann_type) ~ "inframe_ins",
  grepl("splice_acceptor_variant|splice_donor_variant", ann_type) ~ "splice_site",
  grepl("splice_region_variant",  ann_type) ~ "splice_region",
  grepl("exon_loss_variant",      ann_type) ~ "exon_loss",
  grepl("downstream_gene_variant",ann_type) ~ "downstream",
  grepl("upstream_gene_variant",  ann_type) ~ "upstream",
  grepl("intergenic_region",      ann_type) ~ "intergenic",
  grepl("intron_variant",         ann_type) ~ "intron",
  TRUE                                         ~ "other"
)

ann_cat <- factor(
  ann_cat,
  levels = c(
    "stop_gained", "stop_lost", "start_lost",
    "framesshift", "missense", "synonymous",
    "inframe_del", "inframe_ins",
    "splice_site", "splice_region",
    "exon_loss",
    "downstream", "upstream",
    "intron", "intergenic",
    "other"
  )
)

## Little fix: in case of typo in levels
levels(ann_cat)[levels(ann_cat) == "framesshift"] <- "frameshift"

## ============================================================
## 8) SUBSET TO CloneA SAMPLES AND COMPUTE HETEROZYGOSITY
## ============================================================

cloneA_samples <- unique(genomic_types$CloneA)

cat("CloneA samples listed in genomic_type:", length(cloneA_samples), "\n")

samples_in_gds <- samples_all
overlap <- intersect(cloneA_samples, samples_in_gds)

cat("CloneA samples present in GDS:", length(overlap), "\n")
if (length(overlap) == 0) stop("No CloneA sample IDs found in GDS.")

## Filter GDS to CloneA + selected variants
seqSetFilter(
  g,
  variant.id = idx_keep,
  sample.id  = overlap,
  verbose    = FALSE
)

gt <- seqGetData(g, "genotype")   # [ploidy, sample, variant]

het_prop <- apply(gt, 3, function(G) {

  a1 <- G[1, ]
  a2 <- G[2, ]

  nonmiss <- (!is.na(a1)) & (!is.na(a2)) & (a1 >= 0) & (a2 >= 0)

  if (!any(nonmiss, na.rm = TRUE))
    return(NA_real_)

  het <- (a1 != a2) & nonmiss
  sum(het, na.rm = TRUE) / sum(nonmiss, na.rm = TRUE)
})

cat("Heterozygosity summary (CloneA only):\n")
print(summary(het_prop))


het_prop_samp <- apply(gt, 2, function(G) {

  # now G is 2 × n_variants for a single sample
  a1 <- G[1, ]
  a2 <- G[2, ]

  nonmiss <- (!is.na(a1)) & (!is.na(a2)) & (a1 >= 0) & (a2 >= 0)

  if (!any(nonmiss, na.rm = TRUE))
    return(NA_real_)

  het <- (a1 != a2) & nonmiss
  sum(het, na.rm = TRUE) / sum(nonmiss, na.rm = TRUE)
})


## ============================================================
## 9) Build plot_df with region labels (CloneA-based het)
## ============================================================

plot_df <- data.frame(
  chr       = chr_sub,
  pos       = pos_sub,
  ann_type  = ann_type,
  ann_cat   = ann_cat,
  impact    = impact,
  biotype   = biotype,
  het_prop  = het_prop,
  region_id = region_id_sub
)



region_labels <- include_regions %>%
  mutate(region_id = seq_len(n())) %>%
  mutate(label = paste0(chr, ":", win_start, "-", win_end)) %>%
  select(region_id, label)

plot_df <- left_join(plot_df, region_labels, by = "region_id")



het_by_label <- plot_df %>%
  separate(label, into = c("chr", "range"), sep = ":", remove = FALSE) %>%
  separate(range, into = c("start", "end"), sep = "-", convert = TRUE) %>%
  group_by(label, chr, start, end) %>%
  summarise(
    n_variants = n(),
    mean_het   = mean(het_prop, na.rm = TRUE),
    median_het = median(het_prop, na.rm = TRUE),
    sd_het     = sd(het_prop, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(chr, start)

head(het_by_label)

 hetbylabel_plot <- ggplot(het_by_label, aes(x = label, y = mean_het)) +
  geom_col() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))




ggsave(
  file.path(out_dir, "hetbylabel_plot.pdf"),
  plot   = hetbylabel_plot,
  width  = 12,
  height = 8,
  dpi    = 300
)


cat("First rows of plot_df:\n")
print(head(plot_df))

## Some summary stats for CloneA
cat("Effect counts (CloneA variants in regions):\n")
print(sort(table(plot_df$ann_cat), decreasing = TRUE))

cat("Impact counts (CloneA variants in regions):\n")
print(sort(table(plot_df$impact), decreasing = TRUE))

cat("Biotype counts (CloneA variants in regions):\n")
print(sort(table(plot_df$biotype), decreasing = TRUE))

## ============================================================
## 10) Plots + save as PDF (CloneA-based stats)
## ============================================================

## a) Effect types along regions
p_effects <- ggplot(plot_df, aes(x = pos, y = 0, color = ann_cat)) +
  geom_point(alpha = 0.8, position = position_jitter(height = 0.1)) +
  facet_wrap(~ label, scales = "free_x", ncol = 1) +
  theme_bw() +
  theme(
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank()
  ) +
  labs(
    title = "SnpEff effect types across selected regions (CloneA)",
    x     = "Position (bp)",
    color = "Effect category"
  )

ggsave(
  file.path(out_dir, "CloneA_usdobtusa_ann_effects.pdf"),
  plot   = p_effects,
  width  = 12,
  height = 8,
  dpi    = 300
)

## b) Heterozygosity with effect overlay (CloneA)
p_het <- ggplot(plot_df, aes(x = pos, y = het_prop)) +
  geom_line(alpha = 0.4) +
  geom_point(aes(color = ann_cat), size = 1.2) +
  facet_wrap(~ label, scales = "free_x", ncol = 1) +
  theme_bw() +
  labs(
    title = "Per-site heterozygosity with functional annotations (CloneA)",
    x     = "Position (bp)",
    y     = "Proportion heterozygous",
    color = "Effect category"
  )

ggsave(
  file.path(out_dir, "CloneA_usdobtusa_ann_het.pdf"),
  plot   = p_het,
  width  = 12,
  height = 40,
  dpi    = 300
)

## c) Biotypes along regions (CloneA)
p_biotype <- ggplot(plot_df, aes(x = pos, y = 0, color = biotype)) +
  geom_point(alpha = 0.8, position = position_jitter(height = 0.1)) +
  facet_wrap(~ label, scales = "free_x", ncol = 1) +
  theme_bw() +
  theme(
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank()
  ) +
  labs(
    title = "Gene biotypes across selected regions (CloneA)",
    x     = "Position (bp)",
    color = "Biotype"
  )

ggsave(
  file.path(out_dir, "CloneA_usdobtusa_ann_biotypes.pdf"),
  plot   = p_biotype,
  width  = 12,
  height = 8,
  dpi    = 300
)

## d) Effect category counts per region (CloneA)
effect_counts <- plot_df %>%
  count(label, ann_cat)

p_effect_bar <- ggplot(effect_counts, aes(x = label, y = n, fill = ann_cat)) +
  geom_col() +
  coord_flip() +
  theme_bw() +
  labs(
    title = "Effect category counts per region (CloneA)",
    x     = "Region",
    y     = "Variant count",
    fill  = "Effect category"
  )

ggsave(
  file.path(out_dir, "CloneA_usdobtusa_ann_effect_counts_per_region.pdf"),
  plot   = p_effect_bar,
  width  = 12,
  height = 8,
  dpi    = 300
)

## ============================================================
## 11) Export gene list for GO analysis (unfiltered by sample)
## ============================================================

out_genes <- file.path(out_dir, "include_regions_genes_for_GO.txt")

write.table(
  data.frame(gene_id = genes_in_regions),
  file      = out_genes,
  quote     = FALSE,
  sep       = "\t",
  row.names = FALSE,
  col.names = TRUE
)


cat("Wrote gene list for GO to:", out_genes, "\n")

## Optional: long data frame mapping regions and genes (any samples)
gene_region_df <- data.frame(
  chr       = chr_sub,
  pos       = pos_sub,
  region_id = region_id_sub,
  ann       = ann_sub,
  stringsAsFactors = FALSE
)

gene_region_long <- gene_region_df %>%
  mutate(gene_id = gene_list_by_variant) %>%
  unnest(gene_id) %>%
  left_join(region_labels, by = "region_id")

cat("First rows of gene_region_long:\n")
print(head(gene_region_long))

write.csv(gene_region_long, "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/gene_regions.csv")
