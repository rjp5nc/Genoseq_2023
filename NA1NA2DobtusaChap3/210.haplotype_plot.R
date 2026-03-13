library(SeqArray)
library(data.table)
library(ggplot2)

gds_path <- "/scratch/rjp5nc/UK2022_2024/allsites_mito/usdobtusa_mito_allsites.annot2.ALL.gds"

## -----------------------
## 1) Open GDS + basic filtering
## -----------------------
g <- seqOpen(gds_path)

# Keep SNPs only (drop indels) and biallelic only if present
# (Some GDS may not have all of these annotations; we guard with tryCatch.)
try(seqSetFilter(g, variant.sel = seqGetData(g, "variant.id")), silent = TRUE)

# If 'annotation/filter' nodes exist in your GDS, you can add filters here.
# For now, we just pull positions and genotypes and then prune.

pos <- seqGetData(g, "position")
chrom <- seqGetData(g, "chromosome")
vid <- seqGetData(g, "variant.id")
samples <- seqGetData(g, "sample.id")

gt <- seqGetData(g, "genotype")   # array [ploidy, sample, variant]

# haploid mito: use first ploidy slice
gt1 <- gt[1, , ]   # [sample, variant]

# SeqArray genotype coding is typically:
# 0 = missing, 1 = allele1(ref), 2 = allele2(alt), ...
# Convert to: NA (missing), 0 (ref), 1 (alt/anything non-ref)
gt01 <- gt1
gt01[gt01 == 0] <- NA
gt01[!is.na(gt01)] <- as.integer(gt01[!is.na(gt01)] != 1)

## -----------------------
## 2) Keep variable sites (optional but strongly recommended)
## -----------------------
# Drop all-NA sites
keep_non_na <- colSums(!is.na(gt01)) > 0

# Keep sites that vary (not monomorphic among non-missing)
# (At least 2 states observed)
var_sites <- vapply(seq_len(ncol(gt01)), function(j) {
  x <- gt01[, j]
  x <- x[!is.na(x)]
  length(unique(x)) >= 2
}, logical(1))

keep <- keep_non_na & var_sites

gt01 <- gt01[, keep, drop = FALSE]
pos_use <- pos[keep]
vid_use <- vid[keep]
chrom_use <- chrom[keep]

## -----------------------
## 3) Define haplotypes (string across all variable sites)
## -----------------------
# Convert each sample’s 0/1 vector into a string (missing -> N)
hap_str <- apply(gt01, 1, function(x) {
  x2 <- ifelse(is.na(x), "N", as.character(x))
  paste0(x2, collapse = "")
})

hap_id <- as.integer(factor(hap_str))
hap_levels <- data.table(hap_str = unique(hap_str))
hap_counts <- data.table(hap_id = hap_id)[, .N, by = hap_id][order(-N)]

## -----------------------
## 4) Order samples by haplotype similarity (hierarchical clustering)
## -----------------------
# Simple Hamming distance on 0/1 with missing ignored pairwise (quick approximation):
# Convert missing to 0 and compute dist; still good for ordering.
mat_for_dist <- gt01
mat_for_dist[is.na(mat_for_dist)] <- 0
hc <- hclust(dist(mat_for_dist, method = "manhattan"))

sample_order <- samples[hc$order]

## -----------------------
## 5) Long table for plotting
## -----------------------
dt <- as.data.table(gt01)
dt[, sample := samples]
dt[, hap_id := hap_id]
dt_long <- melt(
  dt,
  id.vars = c("sample", "hap_id"),
  variable.name = "var_idx",
  value.name = "allele01"
)

# Map var_idx -> genomic position
dt_long[, var_idx := as.integer(gsub("^V", "", var_idx))]
dt_long[, position := pos_use[var_idx]]

# Order samples for plotting
dt_long[, sample := factor(sample, levels = sample_order)]

# Nice labels for alleles
dt_long[, allele_lab := fifelse(is.na(allele01), "missing",
                         fifelse(allele01 == 0, "ref", "alt"))]

## -----------------------
## 6) Plot: haplotype spanning (tiles across genome)
## -----------------------
p <- ggplot(dt_long, aes(x = position, y = sample, fill = allele_lab)) +
  geom_tile(height = 0.9) +
  scale_fill_manual(values = c(ref = "grey80", alt = "black", missing = "white")) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 5),
    panel.grid = element_blank()
  ) +
  labs(
    title = "Mito haplotype spanning plot (variable sites)",
    x = "Mito position",
    y = "Sample",
    fill = ""
  )

out_png <- "/scratch/rjp5nc/UK2022_2024/allsites_mito/mito_haplotype_spanning.png"
ggsave(out_png, p, width = 60, height = 10, dpi = 300,limitsize = FALSE)

seqClose(g)

message("Saved: ", out_png)
