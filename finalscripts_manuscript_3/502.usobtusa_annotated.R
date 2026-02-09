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
out_dir     <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv"




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


## include_regions should already exist, e.g.:
## include_regions <- data.frame(contig = ..., win_start = ..., win_end = ...)

genomic_types <- read.csv("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/genomic_types_v2.csv")
genomic_types <- subset(genomic_types, Group == "A")

signatures <- read.csv("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/Gilmer_A_het_blocks_signature.csv")

signatures$signature <- sprintf("%08s", signatures$signature)
signatures$signature <- gsub(" ", "0", signatures$signature)

# merged_signatures <- genomic_types %>%
#   left_join(
#     signatures %>%
#       select(sample, signature),
#     by = c("CloneA" = "sample")
#   )

#  genomic_types <- subset(merged_signatures, signature == "01111101")

## ============================================================
## 1) Open GDS and basic variant info
## ============================================================

g <- seqOpen(gds_file)

seqResetFilter(g)

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















library(SeqArray)
library(gdsfmt)
library(SNPRelate)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ape)

## ============================================================
## Paths and region definitions
## ============================================================

gds_file    <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/trimmed10bp_allsites_usdobtusa_renamed_annotated.gds"
contig_file <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/contigmap.txt"
out_dir     <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv"

## Regions to include (original contig names)
include_regions <- data.frame(
  contig = c(
    "JAACYE010000002.1", 
    "JAACYE010000003.1", "JAACYE010000003.1", "JAACYE010000003.1",
    "JAACYE010000007.1", "JAACYE010000007.1",
    "JAACYE010000008.1", "JAACYE010000008.1",
    "JAACYE010000011.1", "JAACYE010000011.1", "JAACYE010000011.1",
    "JAACYE010000012.1"
  ),
  win_start = c(
    4300000,
    2500000, 2600000, 7800000,
    1100000, 1200000, 
    1800000, 1900000, 
    1400000, 1500000, 6500000,
    1100000
  ),
  win_end   = c(
    4399999,
    2599999, 2699999, 7899999,
    1199999, 1299999, 
    1899999, 1999999, 
    1499999, 1599999, 6599999,
    1199999
  )
)

include_regions <- include_regions %>%
  mutate(
    contig    = as.character(contig),
    win_start = as.integer(win_start),
    win_end   = as.integer(win_end),
    region_id = seq_len(n())
  )

## ============================================================
## Genomic types and signatures
## ============================================================

genomic_types <- read.csv("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/genomic_types_v2.csv")
genomic_types <- subset(genomic_types, Group == "A")

signatures <- read.csv("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/Gilmer_A_het_blocks_signature.csv")

## Make signatures 8-digit 0/1 strings
signatures$signature <- sprintf("%08s", signatures$signature)
signatures$signature <- gsub(" ", "0", signatures$signature)

sig_simple <- signatures %>%
  select(sample, signature) %>%
  distinct()

## order of contig_block for mapping bit index (1..nchar(signature)) to blocks
sig_block_order <- signatures %>%
  distinct(contig_block) %>%
  arrange(contig_block) %>%
  pull(contig_block)

cat("Signature block order (for bit indices):\n")
print(sig_block_order)

## CloneA list
cloneA_samples <- unique(genomic_types$CloneA)
cat("CloneA samples listed in genomic_type:", length(cloneA_samples), "\n")

## ============================================================
## Open GDS + basic variant/sample info
## ============================================================

g <- seqOpen(gds_file)
seqResetFilter(g)

contig_all <- seqGetData(g, "chromosome")
pos_all    <- seqGetData(g, "position")
n_var      <- length(pos_all)

samples_all <- seqGetData(g, "sample.id")

cat("Total variants in GDS:", n_var, "\n")
cat("Total samples in GDS :", length(samples_all), "\n")

## Final sample set: CloneA ∩ signatures ∩ GDS
cloneA_sig_samples <- intersect(cloneA_samples, sig_simple$sample)
sample_ids_use     <- intersect(cloneA_sig_samples, samples_all)

cat("CloneA with signatures           :", length(cloneA_sig_samples), "\n")
cat("CloneA with signatures in GDS    :", length(sample_ids_use), "\n")
if (length(sample_ids_use) < 2) stop("Fewer than 2 CloneA+signature samples in GDS.")

## ============================================================
## Map original contigs to GDS contig IDs
## ============================================================

contigmap <- read.table(contig_file, header = FALSE, stringsAsFactors = FALSE)
colnames(contigmap) <- c("old", "new")
map_vec <- setNames(contigmap$new, contigmap$old)

include_regions <- include_regions %>%
  mutate(
    contig_gds = ifelse(
      contig %in% contigmap$old,
      map_vec[contig],
      contig
    )
  )

cat("include_regions (with GDS contig names):\n")
print(include_regions)

## ============================================================
## Define adjacency blocks per contig (merged windows)
## ============================================================

region_blocks <- include_regions %>%
  arrange(contig_gds, win_start) %>%
  group_by(contig_gds, contig) %>%   # keep both contig_gds and original contig
  mutate(
    block_id = cumsum(
      if_else(
        dplyr::row_number() == 1L |
          win_start > dplyr::lag(win_end, default = first(win_end)) + 1L,
        1L, 0L
      )
    )
  ) %>%
  ungroup() %>%
  mutate(
    block_uid    = paste0(contig_gds, "__block", block_id),
    contig_block = sprintf("%s_b%02d", contig, block_id)
  )

cat("# Blocks defined (per contig, merging adjacent/overlapping windows):\n")
print(region_blocks)

## One row per block, with combined start/end coordinates and contig_block ID
block_defs <- region_blocks %>%
  group_by(block_uid, contig_gds, contig, contig_block) %>%
  summarise(
    start = min(win_start),
    end   = max(win_end),
    .groups = "drop"
  ) %>%
  mutate(
    block_label_gds = paste0(contig_gds, ":", start, "-", end),
    block_label_old = paste0(contig,     ":", start, "-", end)
  )

cat("Block definitions:\n")
print(block_defs)

## ============================================================
## Build raw variant index list per block (no SNP filtering yet)
## ============================================================

block_raw_idx_list <- lapply(seq_len(nrow(block_defs)), function(i) {
  ct <- block_defs$contig_gds[i]
  s  <- block_defs$start[i]
  e  <- block_defs$end[i]

  which(contig_all == ct & pos_all >= s & pos_all <= e)
})
names(block_raw_idx_list) <- block_defs$block_uid

cat("Number of raw SNPs per block (before SNPRelate filtering):\n")
print(sapply(block_raw_idx_list, length))

## ============================================================
## Helper: PCA on selected variants (per block), with block-local SNP filtering
## ============================================================



run_region_pca <- function(raw_variant_idx,
                           label,
                           g,
                           sample_ids,
                           maf = 0.00,
                           missing.rate = 0.30) {

  message("Running PCA for ", label, " with ", length(raw_variant_idx), " raw variants")

  if (length(raw_variant_idx) < 2) {
    message("  -> Skipping ", label, ": fewer than 2 raw variants")
    return(NULL)
  }

  ## SNPRelate filtering for THIS block, after selecting the final sample set
  snpset_block <- snpgdsSelectSNP(
    gdsobj        = g,
    sample.id     = sample_ids,
    snp.id        = raw_variant_idx,
    autosome.only = FALSE,
    remove.monosnp = TRUE,
    maf           = maf,
    missing.rate  = missing.rate,
    verbose       = FALSE
  )

  if (length(snpset_block) < 2) {
    message("  -> Skipping ", label, ": <2 SNPs pass MAF/missing filters")
    return(NULL)
  }

  ## Filter GDS to this block + SNP set + samples
  seqSetFilter(
    g,
    variant.id = snpset_block,
    sample.id  = sample_ids,
    verbose    = FALSE
  )

  gt <- seqGetData(g, "genotype")  # [ploidy, sample, variant]

  a1   <- gt[1, , ]
  a2   <- gt[2, , ]
  geno <- a1 + a2

  ## Coerce to matrix if needed
  if (is.null(dim(geno))) {
    geno <- matrix(
      geno,
      nrow = length(sample_ids),
      dimnames = list(sample_ids, NULL)
    )
  }
  rownames(geno) <- sample_ids

  ## Sample filter: drop samples with all missing
  keep_sample <- rowSums(!is.na(geno)) > 0
  if (sum(keep_sample) < 2) {
    message("  -> Skipping ", label, ": fewer than 2 samples with data")
    return(NULL)
  }
  geno <- geno[keep_sample, , drop = FALSE]
  sample_ids_used <- rownames(geno)

  ## Mean impute per SNP
  snp_means <- apply(geno, 2, function(x) mean(x, na.rm = TRUE))
  for (j in seq_len(ncol(geno))) {
    idx_na <- is.na(geno[, j])
    if (any(idx_na)) geno[idx_na, j] <- snp_means[j]
  }

  ## Drop constant SNPs (extra safety)
  snp_sd       <- apply(geno, 2, sd)
  keep_snp_var <- snp_sd > 0 & !is.na(snp_sd)

  if (!any(keep_snp_var)) {
    message("  -> Skipping ", label, ": all SNPs have zero variance after imputation")
    return(NULL)
  }

  geno <- geno[, keep_snp_var, drop = FALSE]

  ## PCA
  pca <- prcomp(geno, center = TRUE, scale. = TRUE)

  scores <- as.data.frame(pca$x)
  scores$sample.id    <- sample_ids_used
  scores$region_label <- label

  var_explained <- (pca$sdev^2) / sum(pca$sdev^2)

  message("  -> PCA OK for ", label, " (", ncol(geno), " SNPs, ",
          nrow(geno), " samples)")

  list(
    scores        = scores,
    var_explained = var_explained
  )
}




i <- 1   # or 2, 3, etc.

label_i <- block_defs$block_label_old[i]
uid_i   <- block_defs$block_uid[i]

raw_idx_i <- block_raw_idx_list[[uid_i]]

seqResetFilter(g)   # always a good idea before running

pca_res_i <- run_region_pca(
  raw_variant_idx = raw_idx_i,
  label           = label_i,
  g               = g,
  sample_ids      = sample_ids_use
)

str(pca_res_i)




## ============================================================
## Helper: NJ tree on selected variants (per block, bit-based coloring)
## ============================================================

run_block_nj <- function(raw_variant_idx,
                         label,
                         contig_block,
                         g,
                         sample_ids,
                         maf = 0.05,
                         missing.rate = 0.10) {

  message("Running NJ tree for ", label, " (", contig_block,
          ") with ", length(raw_variant_idx), " raw variants")

  if (length(raw_variant_idx) < 2) {
    message("  -> Skipping ", label, ": fewer than 2 raw variants")
    return(NULL)
  }

  ## Sample set is already CloneA + signatures
  sample_ids_in <- sample_ids
  if (length(sample_ids_in) < 2) {
    message("  -> Skipping ", label, ": fewer than 2 samples with signatures")
    return(NULL)
  }

  ## SNPRelate filtering for THIS block, after selecting CloneA+signature
  snpset_block <- snpgdsSelectSNP(
    gdsobj        = g,
    sample.id     = sample_ids_in,
    snp.id        = raw_variant_idx,
    autosome.only = FALSE,
    remove.monosnp = TRUE,
    maf           = maf,
    missing.rate  = missing.rate,
    verbose       = FALSE
  )

  if (length(snpset_block) < 2) {
    message("  -> Skipping ", label, ": <2 SNPs pass MAF/missing filters")
    return(NULL)
  }

  ## Filter GDS
  seqSetFilter(
    g,
    variant.id = snpset_block,
    sample.id  = sample_ids_in,
    verbose    = FALSE
  )

  gt <- seqGetData(g, "genotype")

  a1   <- gt[1, , ]
  a2   <- gt[2, , ]
  geno <- a1 + a2

  if (is.null(dim(geno))) {
    geno <- matrix(
      geno,
      nrow = length(sample_ids_in),
      dimnames = list(sample_ids_in, NULL)
    )
  }
  rownames(geno) <- sample_ids_in

  ## Sample filter
  keep_sample <- rowSums(!is.na(geno)) > 0
  if (sum(keep_sample) < 2) {
    message("  -> Skipping ", label, ": fewer than 2 samples with data")
    return(NULL)
  }
  geno <- geno[keep_sample, , drop = FALSE]
  sample_ids_used <- rownames(geno)

  ## Mean impute
  snp_means <- apply(geno, 2, function(x) mean(x, na.rm = TRUE))
  for (j in seq_len(ncol(geno))) {
    idx_na <- is.na(geno[, j])
    if (any(idx_na)) geno[idx_na, j] <- snp_means[j]
  }

  ## Drop constant SNPs
  snp_sd       <- apply(geno, 2, sd)
  keep_snp_var <- snp_sd > 0 & !is.na(snp_sd)

  if (!any(keep_snp_var)) {
    message("  -> Skipping ", label, ": all SNPs have zero variance after imputation")
    return(NULL)
  }

  geno <- geno[, keep_snp_var, drop = FALSE]

  ## Distance + NJ
  d    <- dist(geno)
  tree <- nj(d)
  tree$tip.label <- sample_ids_used

  message("  -> NJ OK for ", label, " (", ncol(geno), " SNPs, ",
          nrow(geno), " samples)")

  list(
    tree         = tree,
    label        = label,
    contig_block = contig_block,
    samples      = sample_ids_used
  )
}

## ============================================================
## Run PCA per block (using CloneA+signature subset)
## ============================================================

seqResetFilter(g)

block_pca_results <- lapply(
  seq_len(nrow(block_defs)),
  function(i) {
    buid  <- block_defs$block_uid[i]
    raw_v <- block_raw_idx_list[[buid]]
    lab   <- block_defs$block_label_old[i]

    run_region_pca(
      raw_variant_idx = raw_v,
      label           = lab,
      g               = g,
      sample_ids      = sample_ids_use
    )
  }
)

block_pca_results <- Filter(Negate(is.null), block_pca_results)

if (length(block_pca_results) > 0) {
  block_pca_scores <- do.call(
    rbind,
    lapply(block_pca_results, `[[`, "scores")
  )

  write.csv(
    block_pca_scores,
    file.path(out_dir, "CloneA_block_PCA_scores_filtered_per_block.csv"),
    row.names = FALSE
  )

  cat("Wrote PCA scores to CloneA_block_PCA_scores_filtered_per_block.csv\n")
} else {
  cat("No block PCA results produced (all blocks skipped).\n")
}






## ============================================================
## Run NJ tree per block (bit-based coloring by contig_block)
## ============================================================

seqResetFilter(g)

block_nj_results <- lapply(
  seq_len(nrow(block_defs)),
  function(i) {
    buid  <- block_defs$block_uid[i]
    raw_v <- block_raw_idx_list[[buid]]
    lab   <- block_defs$block_label_old[i]
    cb    <- block_defs$contig_block[i]

    run_block_nj(
      raw_variant_idx = raw_v,
      label           = lab,
      contig_block    = cb,
      g               = g,
      sample_ids      = sample_ids_use
    )
  }
)

block_nj_results <- Filter(Negate(is.null), block_nj_results)

cat("Number of blocks with NJ trees:", length(block_nj_results), "\n")

## ============================================================
## Plot NJ trees in circular layout, colored by 0/1 at that block's bit
## ============================================================

if (length(block_nj_results) > 0) {

  bit_cols <- c("0" = "dodgerblue3", "1" = "firebrick2")

  n_blocks <- length(block_nj_results)
  n_col    <- 2
  n_row    <- ceiling(n_blocks / n_col)

  pdf(file.path(out_dir, "CloneA_block_NJ_trees_by_block_bit_filtered_per_block.pdf"),
      width = 10, height = 5 * n_row)

  par(mfrow = c(n_row, n_col), mar = c(1, 1, 3, 1))

  for (res in block_nj_results) {
    tree         <- res$tree
    lab          <- res$label
    samples      <- tree$tip.label
    contig_block <- res$contig_block

    ## which bit index corresponds to this block?
    bit_idx <- match(contig_block, sig_block_order)
    if (is.na(bit_idx)) {
      message("contig_block ", contig_block, " not found in sig_block_order; skipping tree for ", lab)
      next
    }

    ## full signatures for these samples
    tip_sig_full <- sig_simple$signature[match(samples, sig_simple$sample)]

    ## extract that bit
    bit_val <- substr(tip_sig_full, bit_idx, bit_idx)

    ## keep only valid 0/1
    keep <- !is.na(bit_val) & bit_val %in% c("0", "1")
    if (sum(keep) < 2) {
      message("  -> Skipping ", lab, ": <2 tips with valid bit values")
      next
    }

    tree_sub    <- drop.tip(tree, which(!keep))
    samples_sub <- samples[keep]
    bit_val     <- bit_val[keep]

    ## reorder bit_val to match tip order after drop.tip
    bit_val <- bit_val[match(tree_sub$tip.label, samples_sub)]

    tip_cols <- bit_cols[bit_val]

    plot.phylo(
      tree_sub,
      type         = "fan",
      tip.color    = tip_cols,
      edge.width   = 0.8,
      label.offset = 0.3,
      cex          = 0.6,
      main         = paste0(lab, "\n", contig_block, " (bit ", bit_idx, ")"),
      font         = 1,
      no.margin    = TRUE
    )
  }

  ## legend for bit 0/1
  plot.new()
  legend(
    "center",
    legend = c("0", "1"),
    col    = bit_cols,
    pch    = 16,
    cex    = 0.9,
    title  = "Signature bit value (per block)"
  )

  dev.off()

  cat("Wrote NJ trees to CloneA_block_NJ_trees_by_block_bit_filtered_per_block.pdf\n")

} else {
  cat("No NJ trees produced (all blocks skipped).\n")
}

## ============================================================
## Clean up
## ============================================================

seqResetFilter(g)


library(dplyr)
library(tidyr)
library(ggplot2)

## ------------------------------------------------------------
## 1) Map each region_label (block) to its contig_block and bit index
## ------------------------------------------------------------

block_map <- block_defs %>%
  dplyr::select(region_label = block_label_old, contig_block) %>%
  mutate(
    bit_idx = match(contig_block, sig_block_order)  # 1..8
  )

## ------------------------------------------------------------
## 2) Attach contig_block, bit_idx, and signature to PCA scores,
##    and compute the bit value (0/1) per sample × block
## ------------------------------------------------------------

block_pca_scores_bit <- block_pca_scores %>%
  left_join(block_map, by = "region_label") %>%
  left_join(sig_simple, by = c("sample.id" = "sample")) %>%
  filter(!is.na(bit_idx), !is.na(signature)) %>%
  mutate(
    bit = substr(signature, bit_idx, bit_idx)
  ) %>%
  filter(bit %in% c("0", "1"))  ## keep only valid 0/1

## If nothing survives for some reason, bail early
if (nrow(block_pca_scores_bit) == 0) {
  stop("No PCA rows with valid bit assignments (check block_map / sig_block_order).")
}

## Order blocks
block_order <- sort(unique(block_pca_scores_bit$region_label))

## ------------------------------------------------------------
## 3) Variance explained per block (PC1–3) from block_pca_results
## ------------------------------------------------------------

var_tab <- lapply(block_pca_results, function(res) {
  lab <- unique(res$scores$region_label)
  data.frame(
    region_label = lab,
    pc           = paste0("PC", seq_along(res$var_explained)),
    pct          = res$var_explained * 100
  )
}) %>%
  bind_rows()

var_wide <- var_tab %>%
  filter(pc %in% c("PC1", "PC2", "PC3")) %>%
  pivot_wider(
    names_from   = pc,
    values_from  = pct,
    names_prefix = "pct_"
  )

## ------------------------------------------------------------
## 4) Build long data frame for plotting
## ------------------------------------------------------------

df_long <- bind_rows(
  block_pca_scores_bit %>%
    select(region_label, bit, PC1, PC2) %>%
    mutate(
      panel = "PC1 vs PC2",
      x     = PC1,
      y     = PC2
    ) %>%
    select(region_label, bit, panel, x, y),

  block_pca_scores_bit %>%
    select(region_label, bit, PC2, PC3) %>%
    mutate(
      panel = "PC2 vs PC3",
      x     = PC2,
      y     = PC3
    ) %>%
    select(region_label, bit, panel, x, y)
) %>%
  left_join(var_wide, by = "region_label") %>%
  mutate(
    region_label = factor(region_label, levels = block_order),
    panel        = factor(panel, levels = c("PC1 vs PC2", "PC2 vs PC3"))
  ) %>%
  filter(!is.na(x), !is.na(y))

## ------------------------------------------------------------
## 5) Annotation df: text with % variance per block & panel
## ------------------------------------------------------------

annotation_df <- var_wide %>%
  mutate(region_label = factor(region_label, levels = block_order)) %>%
  tidyr::expand_grid(
    panel = factor(c("PC1 vs PC2", "PC2 vs PC3"),
                   levels = c("PC1 vs PC2", "PC2 vs PC3"))
  ) %>%
  mutate(
    label = dplyr::case_when(
      panel == "PC1 vs PC2" ~ sprintf("PC1: %.1f%%\nPC2: %.1f%%", pct_PC1, pct_PC2),
      panel == "PC2 vs PC3" ~ sprintf("PC2: %.1f%%\nPC3: %.1f%%", pct_PC2, pct_PC3),
      TRUE ~ ""
    )
  )

## ------------------------------------------------------------
## 6) Plot: rows = blocks, columns = PC1–PC2 / PC2–PC3, colored by bit 0/1
## ------------------------------------------------------------

bit_cols <- c("0" = "dodgerblue3", "1" = "firebrick2")

p_all <- ggplot(df_long, aes(x = x, y = y, color = bit)) +
  geom_point(alpha = 0.7, size = 2) +
  theme_bw() +
  facet_grid(region_label ~ panel, scales = "free") +
  geom_text(
    data        = annotation_df,
    aes(label = label),
    x           = -Inf,
    y           = Inf,
    hjust       = -0.1,
    vjust       = 1.1,
    inherit.aes = FALSE,
    size        = 2.5
  ) +
  scale_color_manual(
    values = bit_cols,
    name   = "Bit value"
  ) +
  labs(
    x     = NULL,
    y     = NULL,
    title = "Per-block PCA (CloneA + signatures): PC1–PC2 and PC2–PC3\nColored by block-specific bit (0/1)"
  ) +
  theme(
    strip.text.y  = element_text(angle = 0),
    panel.spacing = unit(0.5, "lines")
  )

ggsave(
  file.path(out_dir, "PCA_blocks_all_PC1_PC2_PC3_bitcolor.pdf"),
  plot      = p_all,
  width     = 10,
  height    = max(4, 3 * nlevels(df_long$region_label)),
  limitsize = FALSE
)













library(dplyr)
library(stringr)

## ------------------------------------------------------------
## 0) Inputs
## ------------------------------------------------------------

gtf_file <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/Daphnia_obtusa_FS6_genome.gtf"
out_bed  <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/dobtusa_gene_regions.bed"

## ------------------------------------------------------------
## 1) Read GTF
## ------------------------------------------------------------

gtf <- read.delim(
  gtf_file,
  header         = FALSE,
  sep            = "\t",
  comment.char   = "#",
  stringsAsFactors = FALSE
)

colnames(gtf) <- c(
  "seqname", "source", "feature", "start", "end",
  "score", "strand", "frame", "attribute"
)

## Sanity check
# head(gtf, 10)
# table(gtf$feature)

## ------------------------------------------------------------
## 2) Extract gene_id from the attribute column
## ------------------------------------------------------------
## Matches: gene_id DOBTUSG000001;

gtf <- gtf %>%
  mutate(
    gene_id = str_match(attribute, "gene_id\\s+([^;]+);")[, 2]
  )

## Quick check
# head(gtf[, c("feature", "attribute", "gene_id")], 10)

## ------------------------------------------------------------
## 3) Keep only 'gene' features (one row per gene)
## ------------------------------------------------------------

genes <- gtf %>%
  filter(feature == "gene", !is.na(gene_id))

## ------------------------------------------------------------
## 4) Convert to BED format (0-based start, half-open end)
## ------------------------------------------------------------

genes_bed <- genes %>%
  transmute(
    chrom   = seqname,
    start   = as.integer(start) - 1L,  # GTF is 1-based; BED is 0-based
    end     = as.integer(end),
    gene_id = gene_id
  )

## ------------------------------------------------------------
## 5) Write BED file
## ------------------------------------------------------------

write.table(
  genes_bed,
  file      = out_bed,
  quote     = FALSE,
  sep       = "\t",
  row.names = FALSE,
  col.names = FALSE
)

cat("Wrote gene BED to:", out_bed, "\n")














library(SeqArray)
library(dplyr)
library(stringr)
library(GenomicRanges)
library(tibble)

## ------------------------------------------------------------
## 0) Inputs
## ------------------------------------------------------------

gds_file <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/trimmed10bp_allsites_usobtusa_renamed_annotated.gds"
gtf_file <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/Daphnia_obtusa_FS6_genome.gtf"

## ------------------------------------------------------------
## 1) Read GTF and get gene regions (one row per gene)
## ------------------------------------------------------------

gtf <- read.delim(
  gtf_file,
  header       = FALSE,
  sep          = "\t",
  comment.char = "#",
  stringsAsFactors = FALSE
)

colnames(gtf) <- c(
  "seqname", "source", "feature", "start", "end",
  "score", "strand", "frame", "attribute"
)

## Extract gene_id; attributes look like: gene_id DOBTUSG000001; ...
gtf <- gtf %>%
  mutate(
    gene_id = str_match(attribute, "gene_id\\s+([^;]+);")[, 2]
  )

genes <- gtf %>%
  filter(feature == "gene", !is.na(gene_id)) %>%
  mutate(
    start = as.integer(start),
    end   = as.integer(end)
  )

## Build GRanges for gene bodies
gr_genes <- GRanges(
  seqnames = genes$seqname,
  ranges   = IRanges(start = genes$start, end = genes$end),
  gene_id  = genes$gene_id
)

## ------------------------------------------------------------
## 2) Open GDS and get variant positions
## ------------------------------------------------------------

g <- seqOpen(gds_file)
seqResetFilter(g)

chr_all <- seqGetData(g, "chromosome")
pos_all <- seqGetData(g, "position")
var_ids <- seqGetData(g, "variant.id")  # usually 1:n_var
n_var   <- length(pos_all)

cat("Total variants in GDS:", n_var, "\n")

samples_all <- seqGetData(g, "sample.id")
cat("Total samples in GDS:", length(samples_all), "\n")

## (Optional) If you want to reduce memory further, subset samples here:
## seqSetFilter(g, sample.id = some_subset); then refetch samples_all

## ------------------------------------------------------------
## 3) Compute per-site heterozygosity in BLOCKS
## ------------------------------------------------------------

het_site <- rep(NA_real_, n_var)

block_size <- 50000L  # tweak as needed

for (start_idx in seq(1L, n_var, by = block_size)) {
  end_idx <- min(start_idx + block_size - 1L, n_var)

  v_ids_block <- var_ids[start_idx:end_idx]

  ## Set filter to this block of variants (and all samples)
  seqSetFilter(
    g,
    variant.id = v_ids_block,
    sample.id  = samples_all,
    verbose    = FALSE
  )

  gt_block <- seqGetData(g, "genotype")  # [ploidy, sample, block_var]
  # ploidy x n_samples x n_block
  a1 <- gt_block[1, , ]
  a2 <- gt_block[2, , ]

  # Ensure matrices (can drop to vector if only 1 variant)
  if (is.vector(a1)) {
    a1 <- matrix(a1, nrow = length(samples_all))
    a2 <- matrix(a2, nrow = length(samples_all))
  }

  nonmiss <- (!is.na(a1)) & (!is.na(a2)) & (a1 >= 0) & (a2 >= 0)
  het     <- (a1 != a2) & nonmiss

  denom <- colSums(nonmiss, na.rm = TRUE)
  num   <- colSums(het,     na.rm = TRUE)

  het_block <- ifelse(denom > 0, num / denom, NA_real_)

  het_site[start_idx:end_idx] <- het_block

  cat("Processed variants", start_idx, "to", end_idx, "\n")
}

## Clear filters
seqResetFilter(g)

cat("Heterozygosity summary (all sites):\n")
print(summary(het_site))

## ------------------------------------------------------------
## 4) Overlap variants with genes (once)
## ------------------------------------------------------------

gr_snps <- GRanges(
  seqnames = chr_all,
  ranges   = IRanges(start = pos_all, end = pos_all)
)

hits <- findOverlaps(gr_snps, gr_genes)

in_gene <- logical(n_var)
in_gene[queryHits(hits)] <- TRUE
table(in_gene)

stopifnot(length(in_gene) == length(het_site))

## ------------------------------------------------------------
## 5) Genome-wide summaries (optional)
## ------------------------------------------------------------

mean_het_all      <- mean(het_site, na.rm = TRUE)
mean_het_genic    <- mean(het_site[in_gene], na.rm = TRUE)
mean_het_nongenic <- mean(het_site[!in_gene], na.rm = TRUE)

cat("\n=== Mean heterozygosity (site-wise) ===\n")
cat("Whole genome : ", mean_het_all,      "\n")
cat("Genic sites  : ", mean_het_genic,    "\n")
cat("Non-genic    : ", mean_het_nongenic, "\n")

## ------------------------------------------------------------
## 6) Per-gene heterozygosity (no looping over genes)
## ------------------------------------------------------------

hits_df <- as.data.frame(hits)  # queryHits = variant, subjectHits = gene

var_gene_df <- tibble(
  var_idx  = hits_df$queryHits,
  gene_idx = hits_df$subjectHits,
  gene_id  = mcols(gr_genes)$gene_id[gene_idx],
  het      = het_site[var_idx]
)

gene_het <- var_gene_df %>%
  group_by(gene_id) %>%
  summarise(
    n_sites     = n(),
    mean_het    = mean(het, na.rm = TRUE),
    median_het  = median(het, na.rm = TRUE),
    sd_het      = sd(het, na.rm = TRUE),
    min_het     = min(het, na.rm = TRUE),
    max_het     = max(het, na.rm = TRUE),
    .groups     = "drop"
  ) %>%
  left_join(
    genes %>% select(gene_id, seqname, start, end),
    by = "gene_id"
  )

head(gene_het)

## Optionally write out
## write.csv(
##   gene_het,
##   "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/gene_heterozygosity_blockwise.csv",
##   row.names = FALSE
## )

## ------------------------------------------------------------
## 7) Clean up
## ------------------------------------------------------------
