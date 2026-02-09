#ijob -A berglandlab -c10 -p standard --mem=40G
#module load gcc/11.4.0  openmpi/4.1.4 icu R/4.3.1; R

## ============================================================
## NJ tree from mtDNA GDS + two final plots:
##   (1) Full tree, colored by mitotype (rect + circular)
##   (2) Midpoint fan: 1 tip per mitotype-location, tip size = n_ind
##       (drops Unknown mitotypes + UnknownLoc)
## ============================================================

library(SeqArray)
library(SNPRelate)
library(ape)
library(phangorn)
library(ggtree)
library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(tibble)
library(scales)

# ----------------------------
# Paths / outputs
# ----------------------------
gds_path     <- "/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/gatk_gvcf/usdobtusa_mito_joint_Full.gds"
mitotype_csv <- "/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/mito_types.csv"
samples_csv  <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/Sample_ID_Species_merged_20251227.csv"
out_dir      <- "/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/unique_snps_by_group"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# Parameters
# ----------------------------
dp_cutoff   <- 20
size_breaks <- c(1, 5, 25, 60)

# ----------------------------
# Read mitotypes: sample -> mito_plot
# ----------------------------
mito_raw <- read_csv(mitotype_csv, show_col_types = FALSE)
stopifnot(ncol(mito_raw) >= 2)

mito <- mito_raw %>%
  as_tibble() %>%
  { setNames(., c("sample", "mito_type", names(.)[-(1:2)])) } %>%
  transmute(
    sample    = str_trim(as.character(sample)),
    mito_plot = str_trim(as.character(mito_type))
  ) %>%
  filter(nzchar(sample)) %>%
  mutate(mito_plot = ifelse(is.na(mito_plot) | mito_plot == "", "Unknown", mito_plot)) %>%
  distinct(sample, .keep_all = TRUE)

mito_map <- setNames(mito$mito_plot, mito$sample)

# ----------------------------
# Sample list restriction: D. obtusa + NorthAmerica
# ----------------------------
samples <- read.csv(samples_csv)
samples_to_keep <- samples %>%
  filter(Species == "Daphnia obtusa", Continent == "NorthAmerica") %>%
  pull(Sample_ID) %>%
  as.character()

# ----------------------------
# Open GDS + DP filter
# ----------------------------
options(mc.cores = 1)
Sys.setenv(OMP_NUM_THREADS = "1")

g <- seqOpen(gds_path)
on.exit(seqClose(g), add = TRUE)
seqResetFilter(g)

DP <- seqGetData(g, "annotation/format/DP")
if (is.list(DP)) DP <- do.call(cbind, DP)

samps_all <- seqGetData(g, "sample.id")
mean_dp <- if (ncol(DP) == length(samps_all)) colMeans(DP, na.rm = TRUE) else rowMeans(DP, na.rm = TRUE)

keep_samples <- samps_all[!is.na(mean_dp) & mean_dp > dp_cutoff]
keep_samples <- intersect(keep_samples, samples_to_keep)

seqSetFilter(g, sample.id = keep_samples, action = "set")
message("Samples kept after DP + obtusa/NA filter: ", length(keep_samples))

# ----------------------------
# IBS -> NJ tree
# NOTE: snp.id = NULL uses all available SNPs under filter.
# ----------------------------
ibs <- snpgdsIBS(
  g,
  sample.id = keep_samples,
  snp.id = NULL,
  num.thread = 4,
  autosome.only = FALSE
)

D <- 1 - ibs$ibs
dimnames(D) <- list(ibs$sample.id, ibs$sample.id)
tree <- nj(as.dist(D))

# Midpoint root
phy2 <- midpoint(tree)

# Mitotype per tip
tip_mito <- mito_map[phy2$tip.label]
tip_mito[is.na(tip_mito) | tip_mito == ""] <- "Unknown"

# Location mapping:
# Priority: (A) samples Sample_ID -> prefix from Sample_ID_old (SRR -> EBG/PYR/RAP/FS etc)
#          (B) metadata_with_clone Well -> accuratelocation (wins if present)
prefix_only <- function(x) {
  x <- str_trim(as.character(x))
  x <- toupper(x)
  x <- sub("[-_].*$", "", x)
  ifelse(nzchar(x), x, NA_character_)
}

samples_loc_map <- samples %>%
  mutate(
    Sample_ID     = str_trim(as.character(Sample_ID)),
    Sample_ID_old = str_trim(as.character(Sample_ID_old)),
    location      = prefix_only(Sample_ID_old),
    location      = gsub("FS5", "FS", location)
  ) %>%
  filter(!is.na(Sample_ID), Sample_ID != "") %>%
  select(Sample_ID, location) %>%
  distinct(Sample_ID, .keep_all = TRUE)

loc_by_srr <- setNames(samples_loc_map$location, samples_loc_map$Sample_ID)

meta_loc <- metadata_with_clone %>%
  mutate(
    Well = str_trim(as.character(Well)),
    accuratelocation = str_trim(as.character(accuratelocation))
  ) %>%
  filter(!is.na(Well), Well != "", !is.na(accuratelocation), accuratelocation != "")

loc_by_well <- setNames(meta_loc$accuratelocation, meta_loc$Well)

# combine (later wins)
loc_all <- c(loc_by_srr, loc_by_well)
loc_all <- loc_all[!duplicated(names(loc_all), fromLast = TRUE)]

tip_loc <- loc_all[phy2$tip.label]
tip_loc[is.na(tip_loc) | tip_loc == ""] <- "UnknownLoc"

# Drop Unknowns
keep_ok <-  (tip_loc != "UnknownLoc")

#(tip_mito != "Unknown") &


# Build per-tip metadata
tip_meta_all <- tibble(
  tip_index = seq_along(phy2$tip.label),
  label     = phy2$tip.label,
  mito_plot = as.character(tip_mito),
  location  = as.character(tip_loc),
  keep_ok   = keep_ok
) %>%
  mutate(mito_loc = paste0(mito_plot, "_", location))

tip_meta_all <- subset(tip_meta_all, mito_loc != "Unknown_RAP" )
tip_meta_all <- subset(tip_meta_all, mito_loc != "Unknown_BDW3")
tip_meta_all <- subset(tip_meta_all, mito_loc != "Unknown_BDW4")


tip_meta_all$mito_loc <- gsub("Unknown_BDW4", "Unk_BDW", tip_meta_all$mito_loc)

# counts per mito_loc
counts_df <- tip_meta_all %>%
  filter(keep_ok) %>%
  count(mito_loc, mito_plot, location, name = "n_ind")

# representative: first in tree order per mito_loc
rep_tips <- tip_meta_all %>%
  filter(keep_ok) %>%
  group_by(mito_loc) %>%
  slice_min(tip_index, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  left_join(counts_df, by = c("mito_loc", "mito_plot", "location")) %>%
  mutate(label_show = mito_loc)

# keep by unique label to avoid duplicate-label traps
phy2_u <- phy2
phy2_u$tip.label <- paste0(phy2$tip.label, "__i", seq_along(phy2$tip.label))

rep_tips_u <- rep_tips %>%
  mutate(label_u = paste0(label, "__i", tip_index))

phy_rep <- keep.tip(phy2_u, rep_tips_u$label_u)

map_u_to_mitoloc <- setNames(rep_tips_u$mito_loc, rep_tips_u$label_u)
phy_rep$tip.label <- unname(map_u_to_mitoloc[phy_rep$tip.label])

plot_df <- rep_tips_u %>%
  transmute(
    label = mito_loc,
    mito_plot,
    location,
    n_ind,
    label_show
  )

# join sanity
p0 <- ggtree(phy_rep, layout = "fan", open.angle = 0)
stopifnot(setequal(p0$data %>% dplyr::filter(isTip) %>% pull(label),
                   plot_df$label))

p_fan <- p0 %<+% plot_df +
  geom_tippoint(aes(size = n_ind, color = mito_plot), alpha = 0.9) +
  geom_tiplab2(aes(label = label_show, color = mito_plot),
               size = 3.0, offset = 0.02) +
  scale_size_continuous(
    name = "Individuals per mitotype-location",
    breaks = size_breaks
  ) +
  theme_tree() +
  ggtitle("Midpoint-rooted NJ Tree (1 per mitotype-location)")

ggsave(file.path(out_dir, "NJ_midpoint_fan_1perMitotypeLocation_sizedByCount_ggtree.png"),
       p_fan, width = 7, height = 5, dpi = 300)
ggsave(file.path(out_dir, "NJ_midpoint_fan_1perMitotypeLocation_sizedByCount_ggtree.pdf"),
       p_fan, width = 7, height = 5, dpi = 300)


library(patchwork)

load("/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/unique_snps_by_group/superclones_slimmed.RData")

p3 <- p + p_fan


save(
  p_fan, 
  file = file.path(out_dir, "mito_slimmed.RData")
)

ggsave(file.path("/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/unique_snps_by_group/NJ_midpoint_bothmitoandsuper.pdf"),
       p3, width = 18, height = 8, dpi = 300)




























## ============================================================
## Make: allsite_prop_diff_all.png and allsite_prop_diff_all.pdf
## (boxplot of pairwise prop_diff grouped by mitotype comparison)
## ============================================================

library(SeqArray)
library(data.table)
library(readr)
library(stringr)
library(ggplot2)

# ---- Paths ----
gds.fn       <- "/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/usdobtusa_mito_allsites_all.haploid.annotated2.gds"
mitotype_csv <- "/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/mito_types.csv"
samples_csv  <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/Sample_ID_Species_merged_20251227.csv"
out_dir      <- "/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/unique_snps_by_group"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_png <- file.path(out_dir, "allsite_prop_diff_all.png")
out_pdf <- file.path(out_dir, "allsite_prop_diff_all.pdf")

# ---- Parameters (match your pipeline choices) ----
dp_min_mean   <- 30
miss_max_var  <- 0.15
miss_max_samp <- 0.10
keep_outgroup <- "RobertUK_B10"   # included even if it fails filters; set NULL to skip

# ---- Mitotype map: sample.id -> mitotype ----
meta <- fread(mitotype_csv, header = TRUE)
setnames(meta, c("sample.id", "mitotype"))

meta[, sample.id := str_trim(as.character(sample.id))]
meta[, mitotype  := str_trim(as.character(mitotype))]
meta[tolower(mitotype) %in% c("na","nan","none","unknown",""), mitotype := NA_character_]
meta <- meta[!is.na(mitotype) & sample.id != ""]
meta <- meta[, .(mitotype = mitotype[1]), by = sample.id]  # enforce 1:1

# ---- Sample restriction (D. obtusa, remove bdw*) ----
samples <- fread(samples_csv)
samples <- samples[!grepl("^bdw", Sample_ID_old)]
usobtusa_ids <- samples[Species == "Daphnia obtusa", as.character(Sample_ID)]

# ---- Open GDS and apply filters ----
options(mc.cores = 1)
Sys.setenv(OMP_NUM_THREADS = "1")

seqClose(genofile)
g <- seqOpen(gds.fn)
on.exit(seqClose(g), add = TRUE)
seqResetFilter(g)

# Start from D. obtusa sample set
seqSetFilter(g, sample.id = usobtusa_ids, action = "set")

# Missing filters
miss_samp <- seqMissing(g, per.variant = FALSE)
sid_all   <- seqGetData(g, "sample.id")
valid_samp <- sid_all[miss_samp < miss_max_samp]

miss_var <- seqMissing(g, per.variant = TRUE)
vid_all  <- seqGetData(g, "variant.id")
valid_var <- vid_all[miss_var < miss_max_var]

# Mean DP per sample (on current filter)
dp <- seqGetData(g, "annotation/format/DP")
if (is.list(dp)) dp <- do.call(cbind, dp)

# dp can be variants x samples or samples x variants depending on build
# We want mean across variants, per sample:
if (ncol(dp) == length(sid_all)) {
  mean_dp <- colMeans(dp, na.rm = TRUE)
  names(mean_dp) <- sid_all
} else if (nrow(dp) == length(sid_all)) {
  mean_dp <- rowMeans(dp, na.rm = TRUE)
  names(mean_dp) <- sid_all
} else {
  stop("Cannot infer DP orientation: dim(DP) = ", paste(dim(dp), collapse="x"),
       " ; n_samples = ", length(sid_all))
}

dp_keep <- names(mean_dp)[!is.na(mean_dp) & mean_dp > dp_min_mean]

final_samp <- intersect(valid_samp, dp_keep)
final_samp <- intersect(final_samp, usobtusa_ids)

if (!is.null(keep_outgroup)) final_samp <- unique(c(final_samp, keep_outgroup))

seqResetFilter(g)
seqSetFilter(g, sample.id = final_samp, variant.id = valid_var, action = "set")

samps <- seqGetData(g, "sample.id")
message("Final samples: ", length(samps))
message("Final variants: ", length(seqGetData(g, "variant.id")))

# ---- Get genotype matrix (variants x samples), haploid ----
GT_raw <- seqGetData(g, "genotype")
stopifnot(length(dim(GT_raw)) == 3L)

d <- dim(GT_raw)
nsamp <- length(samps)

if (d[2] == nsamp) {
  GT <- t(GT_raw[1, , ])   # variants x samples
} else if (d[3] == nsamp) {
  GT <- GT_raw[1, , ]      # variants x samples
} else {
  stop("Cannot infer genotype dimension order: dim(genotype) = ",
       paste(d, collapse=" x "), " ; nsamp = ", nsamp)
}

GT <- as.matrix(GT)
colnames(GT) <- samps

is_missing <- is.na(GT) | (GT >= 3L)

# ---- Pairwise prop_diff (ALL vs ALL) ----
nsamp <- length(samps)
diff_counts <- matrix(0L, nsamp, nsamp, dimnames = list(samps, samps))
compared    <- matrix(0L, nsamp, nsamp, dimnames = list(samps, samps))

for (i in seq_len(nsamp)) {
  for (j in i:nsamp) {
    ok <- !(is_missing[, i] | is_missing[, j])
    n_ok <- sum(ok)
    compared[i, j] <- compared[j, i] <- n_ok

    if (n_ok > 0) {
      d_ij <- sum(GT[ok, i] != GT[ok, j])
      diff_counts[i, j] <- diff_counts[j, i] <- d_ij
    }
  }
}

diag(diff_counts) <- 0L
diag(compared) <- 0L

prop_diff <- diff_counts / compared

# ---- Long table, deduplicate unordered pairs ----
DT_prop <- as.data.table(as.table(prop_diff))
setnames(DT_prop, c("sampleA","sampleB","prop_diff"))
DT_prop <- DT_prop[sampleA != sampleB]

DT_prop[, sample_min := pmin(sampleA, sampleB)]
DT_prop[, sample_max := pmax(sampleA, sampleB)]
DT_prop <- unique(DT_prop, by = c("sample_min", "sample_max"))

# ---- Annotate mitotypes on both sides ----
pairwise_anno <- merge(
  DT_prop,
  meta[, .(sampleA = sample.id, mitoA = mitotype)],
  by = "sampleA",
  all.x = TRUE
)
pairwise_anno <- merge(
  pairwise_anno,
  meta[, .(sampleB = sample.id, mitoB = mitotype)],
  by = "sampleB",
  all.x = TRUE
)

# drop missing mitotypes
pairwise_anno <- pairwise_anno[!is.na(mitoA) & !is.na(mitoB)]

# undirected mitotype comparison label
pairwise_anno[, mito_comp := paste0(pmin(mitoA, mitoB), "_", pmax(mitoA, mitoB))]


# Define which mitotypes belong to NA1 vs NA2
na1 <- c("D","E","F")
na2 <- c("A","B","C")

pairwise_anno2 <- pairwise_anno %>%
  mutate(
    mitoA = as.character(mitoA),
    mitoB = as.character(mitoB),
    class = case_when(
      mitoA %in% na1 & mitoB %in% na1 ~ "within NA1 (D/E/F)",
      mitoA %in% na2 & mitoB %in% na2 ~ "within NA2 (A/B/C)",
      (mitoA %in% na1 & mitoB %in% na2) | (mitoA %in% na2 & mitoB %in% na1) ~ "between NA1 and NA2",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(class)) %>%
  mutate(class = factor(class,
                        levels = c("within NA2 (A/B/C)", "within NA1 (D/E/F)", "between NA1 and NA2")))

# (optional) keep x labels compact
pairwise_anno2 <- pairwise_anno2 %>%
  mutate(mito_comp = as.character(mito_comp))

p_box2 <- ggplot(pairwise_anno2, aes(x = mito_comp, y = prop_diff)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, height = 0, alpha = 0.6) +
  facet_wrap(~ class, scales = "free_x", nrow = 1) +
  xlab("Mitotype Comparison") +
  ylab("Proportion different (all sites)") +
  ggtitle("Mitotype comparisons (pairwise, all sites)") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "grey90")
  )

ggsave(out_png, p_box2, width = 15, height = 5, dpi = 600)
ggsave(out_pdf, p_box2, width = 15, height = 5, dpi = 600)


allfig <- ( p + p_fan) / p_box2 +
  plot_annotation(tag_levels = "A")

  ggsave("/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/unique_snps_by_group/Trees_and_sim.pdf", allfig, width = 15, height = 10, dpi = 600)
