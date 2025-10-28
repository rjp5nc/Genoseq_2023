#ijob -A berglandlab -c2 -p standard --mem=80G
#module load gcc/11.4.0  openmpi/4.1.4 icu R/4.3.1
#R


library(SNPRelate)
library(SeqArray)
library(gdsfmt)
library(dplyr)
library(ggplot2)
library(tidyr)

# --- Load files ---
genofile.fn <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/trimmed10bp_masked_usobtusa.gds"
gds <- seqOpen(genofile.fn)

metadata <- read.csv("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/2022_2024seqmetadata20250811.csv", header = TRUE)
samplestats <- read.csv("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/usdobtusa_samplestats.csv")
mitotypes <- read.csv("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/mito_types.csv")
genomic_type <- read.csv("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/genomic_types.csv")

# --- Include locations and depth filter ---
include_locations <- c("P66", "P63", "P58", "P62", "Gilmer")

samples_to_keep <- metadata %>%
  inner_join(samplestats, by = c("Well" = "sampleId")) %>%
  filter(accuratelocation %in% include_locations & meanDepth > 3) %>%
  pull(Well)

seqSetFilter(gds, sample.id = samples_to_keep)

# --- SNP filter (example: remove low MAF or high missingness) ---
snpset <- snpgdsSelectSNP(
  gds,
  autosome.only = FALSE,
  remove.monosnp = TRUE,
  maf = 0.05,          # Minimum allele frequency
  missing.rate = 0.1,  # Maximum 10% missingness allowed
  verbose = TRUE
)

# --- Apply filters to the GDS ---
seqSetFilter(gds, sample.id = samples_to_keep, variant.id = snpset)

# --- Identify groups for looping ---
groups <- metadata %>%
  inner_join(genomic_type, by = c("Well" = "CloneA")) %>%
  filter(Well %in% samples_to_keep) %>%
  select(Well, Group)
# --- Prepare output holder ---
pca_all <- data.frame()

# --- Prepare output holders ---
pca_all <- data.frame()
var_summary <- data.frame()

# --- Loop over each group ---
for (g in unique(groups$Group)) {
  cat("Running PCA for group:", g, "\n")
  
  group_samples <- groups %>% filter(Group == g) %>% pull(Well)
  if (length(group_samples) < 3) {
    cat("Skipping group", g, "(too few samples)\n")
    next
  }
  
  seqSetFilter(gds, sample.id = group_samples, variant.id = snpset)
  pca_res <- snpgdsPCA(gds, sample.id = group_samples, num.thread = 10, autosome.only = FALSE)
  
  # Per-sample PCA data
  pca_df <- data.frame(
    sample = pca_res$sample.id,
    PC1 = pca_res$eigenvect[,1],
    PC2 = pca_res$eigenvect[,2],
    PC3 = if (ncol(pca_res$eigenvect) >= 3) pca_res$eigenvect[,3] else rep(0, length(pca_res$sample.id)),
    var_PC1 = pca_res$varprop[1] * 100,
    var_PC2 = pca_res$varprop[2] * 100,
    eigen_PC1 = pca_res$eigenval[1],
    eigen_PC2 = pca_res$eigenval[2],
    Group = g
  )
  
  pca_all <- rbind(pca_all, pca_df)
  
  # Store variance summary per group
  var_summary <- rbind(var_summary, data.frame(
    Group = g,
    eigen_PC1 = pca_res$eigenval[1],
    eigen_PC2 = pca_res$eigenval[2],
    eigen_PC3 = if (length(pca_res$eigenval) >= 3) pca_res$eigenval[3] else NA,
    total_variance = sum(pca_res$eigenval)
  ))
}

# --- Print absolute variances ---
cat("\n--- Total variance per PCA group ---\n")
print(var_summary)

# --- Add facet labels with eigenvalues ---
pca_all <- pca_all %>%
  mutate(
    facet_label = paste0(
      Group, "\n",
      "PC1 var=", signif(eigen_PC1, 3),
      ", PC2 var=", signif(eigen_PC2, 3)
    )
  )

# --- Plot PC1 vs PC2 ---
png("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/pca_usdobtusa_by_group_PC1_PC2_eigenvar.png",
    width = 10000, height = 5000, res = 600)

ggplot(pca_all, aes(x = PC1, y = PC2)) +
  geom_point(size = 2) +
  facet_wrap(~facet_label, scales = "free", ncol = 6) +
  theme_bw() +
  labs(
    title = "Per-Group PCA (Filtered SNPs, Depth > 3)",
    x = "PC1", y = "PC2"
  ) +
  geom_text(
    data = pca_all %>% group_by(Group) %>% slice(1),
    aes(
      x = min(PC1), y = max(PC2),
      label = paste0(
        "Eigenvalues:\nPC1=", signif(eigen_PC1, 3),
        ", PC2=", signif(eigen_PC2, 3)
      )
    ),
    hjust = 0, vjust = -0.5, size = 2.5
  ) +
  theme(strip.text = element_text(size = 8, face = "bold"))

dev.off()






pca_all <- pca_all %>%
  mutate(
    facet_label = paste0(
      Group, "\n",
      "PC1 = ", signif(eigen_PC1, 3), " (", round(var_PC1, 1), "%), ",
      "PC2 = ", signif(eigen_PC2, 3), " (", round(var_PC2, 1), "%)"
    )
  )

# --- Plot PC1 vs PC2 with absolute eigenvalues + % variance ---
png("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/pca_usdobtusa_by_group_PC1_PC2_eigenvar_pct.png",
    width = 10000, height = 5000, res = 600)

ggplot(pca_all, aes(x = PC1, y = PC2)) +
  geom_point(size = 2) +
  facet_wrap(~facet_label, scales = "free", ncol = 6) +
  theme_bw() +
  labs(
    title = "Per-Group PCA (Filtered SNPs, Depth > 3)",
    x = "PC1", y = "PC2"
  ) +
  geom_text(
    data = pca_all %>% group_by(Group) %>% slice(1),
    aes(
      x = min(PC1), y = max(PC2),
      label = paste0(
        "Eigenvalues:\nPC1=", signif(eigen_PC1, 3), " (", round(var_PC1, 1), "%)\n",
        "PC2=", signif(eigen_PC2, 3), " (", round(var_PC2, 1), "%)"
      )
    ),
    hjust = 0, vjust = -0.5, size = 2.5
  ) +
  theme(strip.text = element_text(size = 8, face = "bold"))

dev.off()





















































































































# Full workflow: per-Group PCA -> k-means (k=2..4) on Group PC1-3 -> re-run PCA on each cluster (SNP genotypes)
# Option B behavior: clustering determined from group PCA (PC1-3), but subgroup PCAs are run on SNP genotypes.


library(SNPRelate)
library(SeqArray)
library(gdsfmt)
library(dplyr)
library(ggplot2)
library(tidyr)


# -------------- PARAMETERS ----------------
genofile.fn <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/trimmed10bp_masked_usobtusa.gds"
meta.fn    <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/2022_2024seqmetadata20250811.csv"
samplestats.fn <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/usdobtusa_samplestats.csv"
genomic_type.fn <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/genomic_types.csv"

include_locations <- c("P66","P63","P58","P62","Gilmer")
MIN_DEPTH <- 3
SNP_MAF <- 0.05
SNP_MISS <- 0.1
NUM_THREADS <- 10
SAVE_PLOTS <- TRUE
OUTDIR <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/"
# ------------------------------------------

# load metadata
metadata <- read.csv(meta.fn, header = TRUE, stringsAsFactors = FALSE)
samplestats <- read.csv(samplestats.fn, header = TRUE, stringsAsFactors = FALSE)
genomic_type <- read.csv(genomic_type.fn, header = TRUE, stringsAsFactors = FALSE)

# open GDS using SNPRelate
# gds <- snpgdsOpen(genofile.fn)


# --- Load files ---
gds <- seqOpen(genofile.fn)

# select samples: location + meanDepth > MIN_DEPTH
samples_to_keep <- metadata %>%
  inner_join(samplestats, by = c("Well" = "sampleId")) %>%
  filter(accuratelocation %in% include_locations & meanDepth > MIN_DEPTH) %>%
  pull(Well)

cat("# of selected samples:", length(samples_to_keep), "\n")

# select QC SNPs once (use SNPRelate function)
snp.id <- snpgdsSelectSNP(gds,
                          sample.id = samples_to_keep,
                          autosome.only = FALSE,
                          remove.monosnp = TRUE,
                          maf = SNP_MAF,
                          missing.rate = SNP_MISS,
                          verbose = TRUE)

cat("# of selected SNPs:", length(snp.id), "\n")

# identify groups by joining metadata with genomic_type
groups_tbl <- metadata %>%
  inner_join(genomic_type, by = c("Well" = "CloneA")) %>%
  filter(Well %in% samples_to_keep) %>%
  select(Well, Group)

# containers for summaries
var_summary <- data.frame(Group = character(),
                          eigen_PC1 = numeric(), eigen_PC2 = numeric(), eigen_PC3 = numeric(),
                          total_variance = numeric(),
                          n_samples = integer(),
                          stringsAsFactors = FALSE)

# store per-sample PCA info for plotting/inspection
pca_per_group_samples <- data.frame()

# First: run PCA per group using SNP set (no filtering besides snp.id)
for (g in sort(unique(groups_tbl$Group))) {
  cat("Running group PCA for Group =", g, "\n")
  samples.g <- groups_tbl %>% filter(Group == g) %>% pull(Well)
  n_samps <- length(samples.g)
  if (n_samps < 2) {
    cat(" - skipping (n < 2):", n_samps, "\n")
    next
  }
  # run PCA for the group using SNPRelate; pass the snp.id and sample.id
  pca_g <- snpgdsPCA(gds, sample.id = samples.g, snp.id = snp.id,
                     num.thread = NUM_THREADS, autosome.only = FALSE)
  # eigenvalues and eigenvectors exist
  eigvals <- pca_g$eigenval
  # ensure we have at least 3 eigenvals; pad with NA if missing
  eigen_PC1 <- ifelse(length(eigvals) >= 1, eigvals[1], NA)
  eigen_PC2 <- ifelse(length(eigvals) >= 2, eigvals[2], NA)
  eigen_PC3 <- ifelse(length(eigvals) >= 3, eigvals[3], NA)
  total_var <- sum(eigvals, na.rm = TRUE)
  var_summary <- rbind(var_summary,
                       data.frame(Group = g,
                                  eigen_PC1 = eigen_PC1,
                                  eigen_PC2 = eigen_PC2,
                                  eigen_PC3 = eigen_PC3,
                                  total_variance = total_var,
                                  n_samples = n_samps,
                                  stringsAsFactors = FALSE))
  # store per-sample PC1-3 coordinates for clustering later
  pcs <- pca_g$eigenvect
  # ensure columns exist
  pc1 <- pcs[,1]
  pc2 <- if (ncol(pcs) >= 2) pcs[,2] else rep(0, length(pc1))
  pc3 <- if (ncol(pcs) >= 3) pcs[,3] else rep(0, length(pc1))
  temp <- data.frame(sample = pca_g$sample.id,
                     PC1 = pc1, PC2 = pc2, PC3 = pc3,
                     Group = g, stringsAsFactors = FALSE)
  pca_per_group_samples <- rbind(pca_per_group_samples, temp)

  # optional plot per group (save small figure)
  if (SAVE_PLOTS) {
    pfn <- file.path(OUTDIR, paste0("groupPCA_", g, "_PC1_PC2.png"))
    png(pfn, width = 2000, height = 1600, res = 200)
    plot(temp$PC1, temp$PC2, pch = 19, main = paste0("Group ", g, " PCA (n=", n_samps, ")"),
         xlab = paste0("PC1 (eigen=", signif(eigen_PC1,3), ")"),
         ylab = paste0("PC2 (eigen=", signif(eigen_PC2,3), ")"))
    text(temp$PC1, temp$PC2, labels = temp$sample, cex = 0.6, pos = 3)
    dev.off()
  }
}

# Show groups with total variance
print(var_summary)

subgroup_var_summary <- data.frame(
  Group = character(),
  k = integer(),
  cluster = integer(),
  n_samples = integer(),
  eigen_PC1 = numeric(),
  eigen_PC2 = numeric(),
  eigen_PC3 = numeric(),
  total_variance = numeric(),
  group_total_variance = numeric(),
  stringsAsFactors = FALSE
)

# Filter groups by variance / sample size if desired
var_summary2 <- subset(var_summary, total_variance > 10 & n_samples > 2)
groups_to_split <- var_summary2$Group
cat("Splitting ALL groups:", paste(groups_to_split, collapse=", "), "\n")

for (g in groups_to_split) {
  cat("\nProcessing Group:", g, "\n")
  
  dfg <- pca_per_group_samples %>% filter(Group == g)
  
  # Skip tiny groups
  if (nrow(dfg) < 2) {
    message(sprintf("⚠️ Skipping Group %s: only %d sample(s)", g, nrow(dfg)))
    next
  }
  
  pcs_mat <- as.matrix(dfg %>% select(PC1, PC2, PC3))
  stopifnot(all(is.finite(pcs_mat)))
  
  # Get total variance from original group PCA
  group_total_var <- var_summary %>% filter(Group == g) %>% pull(total_variance)
  
  for (k in 2:6) {
    # Skip impossible k
    if (nrow(dfg) < k) {
      message(sprintf("⚠️ Skipping k=%d for Group %s: n=%d too small", k, g, nrow(dfg)))
      next
    }
    
    cat("  k =", k, " -> running kmeans\n")
    km <- kmeans(pcs_mat, centers = k, nstart = 50)
    dfg$cluster_k <- km$cluster
    
    for (cl in 1:k) {
      samples.cluster <- dfg$sample[dfg$cluster_k == cl]
      n_cl <- length(samples.cluster)
      cat("    cluster", cl, "n=", n_cl, "\n")
      
      # Safe PCA calculation
      if (n_cl >= 2) {
        pca_sub <- tryCatch(
          snpgdsPCA(gds, sample.id = samples.cluster, snp.id = snp.id,
                    num.thread = 1, autosome.only = FALSE),
          error = function(e) return(NULL)
        )
        
        if (!is.null(pca_sub)) {
          eigvals <- pca_sub$eigenval
          eigen1 <- ifelse(length(eigvals) >= 1, eigvals[1], NA)
          eigen2 <- ifelse(length(eigvals) >= 2, eigvals[2], NA)
          eigen3 <- ifelse(length(eigvals) >= 3, eigvals[3], NA)
          total_var <- sum(eigvals, na.rm = TRUE)
        } else {
          eigen1 <- NA; eigen2 <- NA; eigen3 <- NA; total_var <- n_cl
        }
      } else {
        eigen1 <- NA; eigen2 <- NA; eigen3 <- NA; total_var <- n_cl
      }
      
      subgroup_var_summary <- rbind(subgroup_var_summary,
                                    data.frame(
                                      Group = g,
                                      k = k,
                                      cluster = cl,
                                      n_samples = n_cl,
                                      eigen_PC1 = eigen1,
                                      eigen_PC2 = eigen2,
                                      eigen_PC3 = eigen3,
                                      total_variance = total_var,
                                      group_total_variance = group_total_var
                                    ))
    }
  }
}

# Save results
write.csv(var_summary,
          file = file.path(OUTDIR, "group_eigen_summary.csv"),
          row.names = FALSE)
write.csv(subgroup_var_summary,
          file = file.path(OUTDIR, "subgroup_eigen_summary.csv"),
          row.names = FALSE)

cat("\n✅ FINISHED! Subgroup variance summary:\n")
print(subgroup_var_summary)


snpgdsClose(gds)


























library(SNPRelate)
library(SeqArray)
library(gdsfmt)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggrepel)


# -------------- PARAMETERS ----------------
genofile.fn <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/trimmed10bp_masked_usobtusa.gds"
meta.fn    <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/2022_2024seqmetadata20250811.csv"
samplestats.fn <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/usdobtusa_samplestats.csv"
genomic_type.fn <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/genomic_types.csv"

include_locations <- c("P66","P63","P58","P62","Gilmer")
MIN_DEPTH <- 3
SNP_MAF <- 0.05
SNP_MISS <- 0.1
NUM_THREADS <- 10
SAVE_PLOTS <- TRUE
OUTDIR <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/"
# ------------------------------------------

# load metadata
metadata <- read.csv(meta.fn, header = TRUE, stringsAsFactors = FALSE)
samplestats <- read.csv(samplestats.fn, header = TRUE, stringsAsFactors = FALSE)
genomic_type <- read.csv(genomic_type.fn, header = TRUE, stringsAsFactors = FALSE)

# open GDS using SNPRelate
# gds <- snpgdsOpen(genofile.fn)


# --- Load files ---
gds <- seqOpen(genofile.fn)

# select samples: location + meanDepth > MIN_DEPTH
samples_to_keep <- metadata %>%
  inner_join(samplestats, by = c("Well" = "sampleId")) %>%
  filter(accuratelocation %in% include_locations & meanDepth > MIN_DEPTH) %>%
  pull(Well)

cat("# of selected samples:", length(samples_to_keep), "\n")

# select QC SNPs once (use SNPRelate function)
snp.id <- snpgdsSelectSNP(gds,
                          sample.id = samples_to_keep,
                          autosome.only = FALSE,
                          remove.monosnp = TRUE,
                          maf = SNP_MAF,
                          missing.rate = SNP_MISS,
                          verbose = TRUE)

cat("# of selected SNPs:", length(snp.id), "\n")

# identify groups by joining metadata with genomic_type
groups_tbl <- metadata %>%
  inner_join(genomic_type, by = c("Well" = "CloneA")) %>%
  filter(Well %in% samples_to_keep) %>%
  select(Well, Group)

# containers for summaries
var_summary <- data.frame(Group = character(),
                          eigen_PC1 = numeric(), eigen_PC2 = numeric(), eigen_PC3 = numeric(),
                          total_variance = numeric(),
                          n_samples = integer(),
                          stringsAsFactors = FALSE)

# store per-sample PCA info for plotting/inspection
pca_per_group_samples <- data.frame()

# First: run PCA per group using SNP set (no filtering besides snp.id)
for (g in sort(unique(groups_tbl$Group))) {
  cat("Running group PCA for Group =", g, "\n")
  samples.g <- groups_tbl %>% filter(Group == g) %>% pull(Well)
  n_samps <- length(samples.g)
  if (n_samps < 2) {
    cat(" - skipping (n < 2):", n_samps, "\n")
    next
  }
  # run PCA for the group using SNPRelate; pass the snp.id and sample.id
  pca_g <- snpgdsPCA(gds, sample.id = samples.g, snp.id = snp.id,
                     num.thread = NUM_THREADS, autosome.only = FALSE)
  # eigenvalues and eigenvectors exist
  eigvals <- pca_g$eigenval
  # ensure we have at least 3 eigenvals; pad with NA if missing
  eigen_PC1 <- ifelse(length(eigvals) >= 1, eigvals[1], NA)
  eigen_PC2 <- ifelse(length(eigvals) >= 2, eigvals[2], NA)
  eigen_PC3 <- ifelse(length(eigvals) >= 3, eigvals[3], NA)
  total_var <- sum(eigvals, na.rm = TRUE)
  var_summary <- rbind(var_summary,
                       data.frame(Group = g,
                                  eigen_PC1 = eigen_PC1,
                                  eigen_PC2 = eigen_PC2,
                                  eigen_PC3 = eigen_PC3,
                                  total_variance = total_var,
                                  n_samples = n_samps,
                                  stringsAsFactors = FALSE))
  # store per-sample PC1-3 coordinates for clustering later
  pcs <- pca_g$eigenvect
  # ensure columns exist
  pc1 <- pcs[,1]
  pc2 <- if (ncol(pcs) >= 2) pcs[,2] else rep(0, length(pc1))
  pc3 <- if (ncol(pcs) >= 3) pcs[,3] else rep(0, length(pc1))
  temp <- data.frame(sample = pca_g$sample.id,
                     PC1 = pc1, PC2 = pc2, PC3 = pc3,
                     Group = g, stringsAsFactors = FALSE)
  pca_per_group_samples <- rbind(pca_per_group_samples, temp)

  # optional plot per group (save small figure)
  if (SAVE_PLOTS) {
    pfn <- file.path(OUTDIR, paste0("groupPCA_", g, "_PC1_PC2.png"))
    png(pfn, width = 2000, height = 1600, res = 200)
    plot(temp$PC1, temp$PC2, pch = 19, main = paste0("Group ", g, " PCA (n=", n_samps, ")"),
         xlab = paste0("PC1 (eigen=", signif(eigen_PC1,3), ")"),
         ylab = paste0("PC2 (eigen=", signif(eigen_PC2,3), ")"))
    text(temp$PC1, temp$PC2, labels = temp$sample, cex = 0.6, pos = 3)
    dev.off()
  }
}

# Show groups with total variance
print(var_summary)

cluster_assignments <- data.frame(
  sample = character(), Group = character(), k = integer(), cluster = integer(),
  PC1 = numeric(), PC2 = numeric(), PC3 = numeric(), stringsAsFactors = FALSE
)

# We'll keep your previous subgroup_var_summary population, but now also record per-sample clusters
subgroup_var_summary <- data.frame(
  Group = character(), k = integer(), cluster = integer(), n_samples = integer(),
  eigen_PC1 = numeric(), eigen_PC2 = numeric(), eigen_PC3 = numeric(),
  total_variance = numeric(), group_total_variance = numeric(), stringsAsFactors = FALSE
)

# Filter groups by variance / sample size if desired (as before)
var_summary2 <- subset(var_summary, total_variance > 10 & n_samples > 2)
groups_to_split <- var_summary2$Group
cat("Splitting ALL groups:", paste(groups_to_split, collapse=", "), "\n")



for (g in groups_to_split) {
  cat("\nProcessing Group:", g, "\n")
  dfg_base <- pca_per_group_samples %>% filter(Group == g)

  if (nrow(dfg_base) < 2) {
    message(sprintf("⚠️ Skipping Group %s: only %d sample(s)", g, nrow(dfg_base)))
    next
  }

  pcs_mat <- as.matrix(dfg_base %>% select(PC1, PC2, PC3))
  stopifnot(all(is.finite(pcs_mat)))

  group_total_var <- var_summary %>% filter(Group == g) %>% pull(total_variance)

  for (k in 2:6) {  # ~5 groupings: k = 2,3,4,5,6
    if (nrow(dfg_base) < k) {
      message(sprintf("⚠️ Skipping k=%d for Group %s: n=%d too small", k, g, nrow(dfg_base)))
      next
    }

    cat("  k =", k, " -> running kmeans\n")
    km <- kmeans(pcs_mat, centers = k, nstart = 50)

    dfg_k <- dfg_base %>%
      mutate(cluster = km$cluster, k = k) %>%
      select(sample, Group, k, cluster, PC1, PC2, PC3)

    # store per-sample assignments for plotting
    cluster_assignments <- bind_rows(cluster_assignments, dfg_k)

    # summarize PCA-on-genotypes per cluster (as you had)
    for (cl in 1:k) {
      samples.cluster <- dfg_k$sample[dfg_k$cluster == cl]
      n_cl <- length(samples.cluster)
      cat("    cluster", cl, "n=", n_cl, "\n")

      if (n_cl >= 2) {
        pca_sub <- tryCatch(
          snpgdsPCA(gds, sample.id = samples.cluster, snp.id = snp.id,
                    num.thread = 1, autosome.only = FALSE),
          error = function(e) NULL
        )

        if (!is.null(pca_sub)) {
          eigvals <- pca_sub$eigenval
          eigen1 <- ifelse(length(eigvals) >= 1, eigvals[1], NA)
          eigen2 <- ifelse(length(eigvals) >= 2, eigvals[2], NA)
          eigen3 <- ifelse(length(eigvals) >= 3, eigvals[3], NA)
          total_var <- sum(eigvals, na.rm = TRUE)
        } else {
          eigen1 <- NA; eigen2 <- NA; eigen3 <- NA; total_var <- n_cl
        }
      } else {
        eigen1 <- NA; eigen2 <- NA; eigen3 <- NA; total_var <- n_cl
      }

      subgroup_var_summary <- rbind(subgroup_var_summary, data.frame(
        Group = g, k = k, cluster = cl, n_samples = n_cl,
        eigen_PC1 = eigen1, eigen_PC2 = eigen2, eigen_PC3 = eigen3,
        total_variance = total_var, group_total_variance = group_total_var
      ))
    }
  }
}

# ----------------- ADD: make per-Group multi-panel figures (PC1 vs PC2) -----------------
if (SAVE_PLOTS && nrow(cluster_assignments) > 0) {
  for (g in sort(unique(cluster_assignments$Group))) {
    dfg <- cluster_assignments %>% filter(Group == g)

    # only keep k values that actually ran for this group
    dfg$k <- factor(dfg$k, levels = sort(unique(dfg$k)))

    # Use fixed scales so panels are comparable within a group
    p <- ggplot(dfg, aes(x = PC1, y = PC2, color = as.factor(cluster))) +
      geom_point(size = 2, alpha = 0.85) +
      ggrepel::geom_text_repel(aes(label = sample), size = 2.4, max.overlaps = 60, show.legend = FALSE) +
      facet_wrap(~ k, ncol = 3, scales = "fixed") +
      labs(
        title = paste0("Group ", g, " — PC1 vs PC2 with k-means clusters"),
        subtitle = "Clustering performed on Group PC1–PC3; points are Group PCA scores (PC1, PC2)",
        x = "PC1 (Group PCA)", y = "PC2 (Group PCA)", color = "Cluster"
      ) +
      theme_bw(base_size = 12) +
      theme(
        legend.position = "bottom",
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold")
      )
    OUTDIR <- sub("/+$", "", OUTDIR)
    pfn <- file.path(OUTDIR,  paste0("Group_", g, "_PC1PC2_k2to6.png"))
    ggsave(pfn, plot = p, width = 12, height = 8, dpi = 300)
    message("Saved: ", pfn)
  }
}

# ----------------- keep your CSV saves (if not already run later) -----------------
write.csv(var_summary, file = file.path(OUTDIR, "group_eigen_summary.csv"), row.names = FALSE)
write.csv(subgroup_var_summary, file = file.path(OUTDIR, "subgroup_eigen_summary.csv"), row.names = FALSE)






































library(SNPRelate)
library(SeqArray)
library(gdsfmt)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggrepel)


# -------------- PARAMETERS ----------------
genofile.fn <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/trimmed10bp_masked_usobtusa.gds"
meta.fn    <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/2022_2024seqmetadata20250811.csv"
samplestats.fn <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/usdobtusa_samplestats.csv"
genomic_type.fn <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/genomic_types.csv"

include_locations <- c("P66","P63","P58","P62","Gilmer")
MIN_DEPTH <- 3
SNP_MAF <- 0.05
SNP_MISS <- 0.1
NUM_THREADS <- 10
SAVE_PLOTS <- TRUE
OUTDIR <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/"
# ------------------------------------------

# load metadata
metadata <- read.csv(meta.fn, header = TRUE, stringsAsFactors = FALSE)
samplestats <- read.csv(samplestats.fn, header = TRUE, stringsAsFactors = FALSE)
genomic_type <- read.csv(genomic_type.fn, header = TRUE, stringsAsFactors = FALSE)

# open GDS using SNPRelate
# gds <- snpgdsOpen(genofile.fn)

# --- Load files ---
gds <- seqOpen(genofile.fn)

# select samples: location + meanDepth > MIN_DEPTH
samples_to_keep <- metadata %>%
  inner_join(samplestats, by = c("Well" = "sampleId")) %>%
  filter(accuratelocation %in% include_locations & meanDepth > MIN_DEPTH) %>%
  pull(Well)

cat("# of selected samples:", length(samples_to_keep), "\n")

# select QC SNPs once (use SNPRelate function)
snp.id <- snpgdsSelectSNP(gds,
                          sample.id = samples_to_keep,
                          autosome.only = FALSE,
                          remove.monosnp = TRUE,
                          maf = SNP_MAF,
                          missing.rate = SNP_MISS,
                          verbose = TRUE)

cat("# of selected SNPs:", length(snp.id), "\n")

# identify groups by joining metadata with genomic_type
groups_tbl <- metadata %>%
  inner_join(genomic_type, by = c("Well" = "CloneA")) %>%
  filter(Well %in% samples_to_keep) %>%
  select(Well, Group)

# containers for summaries
var_summary <- data.frame(Group = character(),
                          eigen_PC1 = numeric(), eigen_PC2 = numeric(), eigen_PC3 = numeric(),
                          total_variance = numeric(),
                          n_samples = integer(),
                          stringsAsFactors = FALSE)

# store per-sample PCA info for plotting/inspection
pca_per_group_samples <- data.frame()

# First: run PCA per group using SNP set (no filtering besides snp.id)
for (g in sort(unique(groups_tbl$Group))) {
  cat("Running group PCA for Group =", g, "\n")
  samples.g <- groups_tbl %>% filter(Group == g) %>% pull(Well)
  n_samps <- length(samples.g)
  if (n_samps < 2) {
    cat(" - skipping (n < 2):", n_samps, "\n")
    next
  }
  # run PCA for the group using SNPRelate; pass the snp.id and sample.id
  pca_g <- snpgdsPCA(gds, sample.id = samples.g, snp.id = snp.id,
                     num.thread = NUM_THREADS, autosome.only = FALSE)
  # eigenvalues and eigenvectors exist
  eigvals <- pca_g$eigenval
  # ensure we have at least 3 eigenvals; pad with NA if missing
  eigen_PC1 <- ifelse(length(eigvals) >= 1, eigvals[1], NA)
  eigen_PC2 <- ifelse(length(eigvals) >= 2, eigvals[2], NA)
  eigen_PC3 <- ifelse(length(eigvals) >= 3, eigvals[3], NA)
  total_var <- sum(eigvals, na.rm = TRUE)
  var_summary <- rbind(var_summary,
                       data.frame(Group = g,
                                  eigen_PC1 = eigen_PC1,
                                  eigen_PC2 = eigen_PC2,
                                  eigen_PC3 = eigen_PC3,
                                  total_variance = total_var,
                                  n_samples = n_samps,
                                  stringsAsFactors = FALSE))
  # store per-sample PC1-3 coordinates for clustering later
  pcs <- pca_g$eigenvect
  # ensure columns exist
  pc1 <- pcs[,1]
  pc2 <- if (ncol(pcs) >= 2) pcs[,2] else rep(0, length(pc1))
  pc3 <- if (ncol(pcs) >= 3) pcs[,3] else rep(0, length(pc1))
  temp <- data.frame(sample = pca_g$sample.id,
                     PC1 = pc1, PC2 = pc2, PC3 = pc3,
                     Group = g, stringsAsFactors = FALSE)
  pca_per_group_samples <- rbind(pca_per_group_samples, temp)

  # optional plot per group (save small figure)
  if (SAVE_PLOTS) {
    pfn <- file.path(OUTDIR, paste0("groupPCA_", g, "_PC1_PC2.png"))
    png(pfn, width = 2000, height = 1600, res = 200)
    plot(temp$PC1, temp$PC2, pch = 19, main = paste0("Group ", g, " PCA (n=", n_samps, ")"),
         xlab = paste0("PC1 (eigen=", signif(eigen_PC1,3), ")"),
         ylab = paste0("PC2 (eigen=", signif(eigen_PC2,3), ")"))
    text(temp$PC1, temp$PC2, labels = temp$sample, cex = 0.6, pos = 3)
    dev.off()
  }
}

OUTDIR <- sub("/+$", "", OUTDIR)
if (!dir.exists(OUTDIR)) dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

subgroup_var_summary <- data.frame(
  Group = character(),
  k = integer(),
  cluster = integer(),
  n_samples = integer(),
  eigen_PC1 = numeric(),
  eigen_PC2 = numeric(),
  eigen_PC3 = numeric(),
  total_variance = numeric(),
  group_total_variance = numeric(),
  stringsAsFactors = FALSE
)

var_summary2 <- subset(var_summary, total_variance > 10 & n_samples > 2)
groups_to_split <- var_summary2$Group
cat("Splitting ALL groups:", paste(groups_to_split, collapse = ", "), "\n")

all_plot_dfs <- list()  # container before loop


for (g in groups_to_split) {
  cat("\nProcessing Group:", g, "\n")

  dfg_base <- pca_per_group_samples %>% filter(Group == g)
  if (nrow(dfg_base) < 2) {
    message(sprintf("⚠️ Skipping Group %s: only %d sample(s)", g, nrow(dfg_base)))
    next
  }

  pcs_mat <- as.matrix(dfg_base %>% select(PC1, PC2, PC3))
  stopifnot(all(is.finite(pcs_mat)))
  group_total_var <- var_summary %>% filter(Group == g) %>% pull(total_variance)

  # collect per-k assignments for THIS group only
  plot_df <- dfg_base[0, c("sample","Group","PC1","PC2","PC3")]
  plot_df$k <- integer(0)
  plot_df$cluster <- integer(0)

  for (k in 2:6) {
    if (nrow(dfg_base) < k) {
      message(sprintf("⚠️ Skipping k=%d for Group %s: n=%d too small", k, g, nrow(dfg_base)))
      next
    }

    cat("  k =", k, " -> running kmeans\n")
    km <- kmeans(pcs_mat, centers = k, nstart = 50)
    dfg_k <- dfg_base %>%
      mutate(k = k, cluster = km$cluster) %>%
      select(sample, Group, k, cluster, PC1, PC2, PC3)

    # append to this group's plotting df
    plot_df <- bind_rows(plot_df, dfg_k)

    # your PCA-on-genotypes summaries per cluster
    for (cl in 1:k) {
      samples.cluster <- dfg_k$sample[dfg_k$cluster == cl]
      n_cl <- length(samples.cluster)
      cat("    cluster", cl, "n=", n_cl, "\n")

      if (n_cl >= 2) {
        pca_sub <- tryCatch(
          snpgdsPCA(gds, sample.id = samples.cluster, snp.id = snp.id,
                    num.thread = 1, autosome.only = FALSE),
          error = function(e) NULL
        )
        if (!is.null(pca_sub)) {
          eigvals <- pca_sub$eigenval
          eigen1 <- ifelse(length(eigvals) >= 1, eigvals[1], NA)
          eigen2 <- ifelse(length(eigvals) >= 2, eigvals[2], NA)
          eigen3 <- ifelse(length(eigvals) >= 3, eigvals[3], NA)
          total_var <- sum(eigvals, na.rm = TRUE)
        } else {
          eigen1 <- NA; eigen2 <- NA; eigen3 <- NA; total_var <- n_cl
        }
      } else {
        eigen1 <- NA; eigen2 <- NA; eigen3 <- NA; total_var <- n_cl
      }

      subgroup_var_summary <- rbind(subgroup_var_summary, data.frame(
        Group = g, k = k, cluster = cl, n_samples = n_cl,
        eigen_PC1 = eigen1, eigen_PC2 = eigen2, eigen_PC3 = eigen3,
        total_variance = total_var, group_total_variance = group_total_var,
        stringsAsFactors = FALSE
      ))
    }
  } # end k loop

  # ---- MAKE & SAVE PLOT FOR THIS GROUP RIGHT HERE ----
  plot_group_clusters <- function(g) {
  dfg <- subset(pca_per_group_samples, Group == g)
  n <- nrow(dfg)
  if (n < 2) { message(sprintf("Skipping %s: only %d samples", g, n)); return(invisible(NULL)) }

  pcs <- as.matrix(dfg[, c("PC1","PC2","PC3")])
  stopifnot(all(is.finite(pcs)))

  # Deduplicate PC rows to avoid k > distinct rows
  key <- apply(round(pcs, 10), 1, paste, collapse = "|")
  uniq_idx   <- !duplicated(key)
  Xuniq      <- pcs[uniq_idx, , drop = FALSE]
  uniq_keys  <- key[uniq_idx]
  distinct_n <- nrow(Xuniq)

  # Candidate ks: 2..6, but must be <= sample count AND <= distinct PC rows
  k_candidates <- intersect(2:6, 2:min(n, distinct_n))
  if (length(k_candidates) == 0) {
    message(sprintf("[Group %s] no valid k (n=%d, distinct=%d) — skipping", g, n, distinct_n))
    return(invisible(NULL))
  }

  rows <- list()
  for (k in k_candidates) {
    if (k > nrow(Xuniq)) {
      message(sprintf("[Group %s] skip k=%d: only %d distinct rows", g, k, nrow(Xuniq)))
      next
    }

    km <- tryCatch(kmeans(Xuniq, centers = k, nstart = 50),
                   error = function(e) { 
                     message(sprintf("[Group %s] k=%d failed: %s", g, k, e$message)); 
                     return(NULL)
                   })
    if (is.null(km)) next

    cl_uniq <- km$cluster
    # Ensure no empty clusters
    if (any(table(cl_uniq) == 0L) || length(unique(cl_uniq)) < 2L) {
      message(sprintf("[Group %s] skip k=%d: collapsed/empty clusters", g, k))
      next
    }

    # Map unique-row cluster back to all samples (duplicates get same label)
    cl_all <- cl_uniq[match(key, uniq_keys)]

    rows[[as.character(k)]] <- transform(
      dfg[, c("sample","Group","PC1","PC2","PC3")],
      k = k, cluster = cl_all
    )
  }

  if (length(rows) == 0) {
    message(sprintf("[Group %s] nothing to plot after filtering invalid ks", g))
    return(invisible(NULL))
  }

  plot_df <- dplyr::bind_rows(rows)
  plot_df$k <- factor(plot_df$k, levels = sort(unique(plot_df$k)))

  use_labels <- nrow(plot_df) <= 200
  p <- ggplot(plot_df, aes(PC1, PC2, color = as.factor(cluster))) +
    geom_point(size = 2, alpha = 0.9) +
    facet_wrap(~ k, ncol = 3, scales = "fixed") +
    coord_equal() +
    labs(
      title = paste0("Group ", g, " - PC1 vs PC2 with k-means clusters"),
      subtitle = "Clustering on unique rows of Group PC1–PC3; points are Group PCA PC1–PC2",
      x = "PC1 (Group PCA)", y = "PC2 (Group PCA)", color = "Cluster"
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom",
          panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"))

  pfn <- file.path(OUTDIR, paste0("Group_", g, "_PC1PC2_k2to6.png"))
  message("Saving -> ", pfn); ggsave(pfn, plot = p, width = 12, height = 8, dpi = 300)
  message("Saved: ", pfn)
  invisible(NULL)


   if (nrow(plot_df) > 0) {
    # ensure stable column order and types
    plot_df <- plot_df %>%
      mutate(Group = as.character(Group),
             k = as.integer(k),
             cluster = as.integer(cluster)) %>%
      select(Group, sample, k, cluster, PC1, PC2, PC3)
    all_plot_dfs[[as.character(g)]] <- plot_df
  } else {
    message(sprintf("[Group %s] no clustering rows to collect", g))
  }
  }
}

















































library(SNPRelate)
library(SeqArray)
library(gdsfmt)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggrepel)

# -------------- PARAMETERS ----------------
genofile.fn <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/trimmed10bp_masked_usobtusa.gds"
meta.fn    <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/2022_2024seqmetadata20250811.csv"
samplestats.fn <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/usdobtusa_samplestats.csv"
genomic_type.fn <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/genomic_types.csv"

include_locations <- c("P66","P63","P58","P62","Gilmer")
MIN_DEPTH <- 3
SNP_MAF <- 0.05
SNP_MISS <- 0.1
NUM_THREADS <- 10
SAVE_PLOTS <- TRUE
OUTDIR <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/"
# ------------------------------------------

# load metadata
metadata <- read.csv(meta.fn, header = TRUE, stringsAsFactors = FALSE)
samplestats <- read.csv(samplestats.fn, header = TRUE, stringsAsFactors = FALSE)
genomic_type <- read.csv(genomic_type.fn, header = TRUE, stringsAsFactors = FALSE)

# --- Load files ---
gds <- seqOpen(genofile.fn)

# select samples: location + meanDepth > MIN_DEPTH
samples_to_keep <- metadata %>%
  inner_join(samplestats, by = c("Well" = "sampleId")) %>%
  filter(accuratelocation %in% include_locations & meanDepth > MIN_DEPTH) %>%
  pull(Well)

cat("# of selected samples:", length(samples_to_keep), "\n")

# select QC SNPs once (use SNPRelate function)
snp.id <- snpgdsSelectSNP(gds,
                          sample.id = samples_to_keep,
                          autosome.only = FALSE,
                          remove.monosnp = TRUE,
                          maf = SNP_MAF,
                          missing.rate = SNP_MISS,
                          verbose = TRUE)

cat("# of selected SNPs:", length(snp.id), "\n")

# identify groups by joining metadata with genomic_type
groups_tbl <- metadata %>%
  inner_join(genomic_type, by = c("Well" = "CloneA")) %>%
  filter(Well %in% samples_to_keep) %>%
  select(Well, Group)

# containers for summaries
var_summary <- data.frame(Group = character(),
                          eigen_PC1 = numeric(), eigen_PC2 = numeric(), eigen_PC3 = numeric(),
                          total_variance = numeric(),
                          n_samples = integer(),
                          stringsAsFactors = FALSE)

# store per-sample PCA info for plotting/inspection
pca_per_group_samples <- data.frame()

# First: run PCA per group using SNP set (no filtering besides snp.id)
for (g in sort(unique(groups_tbl$Group))) {
  cat("Running group PCA for Group =", g, "\n")
  samples.g <- groups_tbl %>% filter(Group == g) %>% pull(Well)
  n_samps <- length(samples.g)
  if (n_samps < 2) {
    cat(" - skipping (n < 2):", n_samps, "\n")
    next
  }
  # run PCA for the group using SNPRelate; pass the snp.id and sample.id
  pca_g <- snpgdsPCA(gds, sample.id = samples.g, snp.id = snp.id,
                     num.thread = NUM_THREADS, autosome.only = FALSE)
  # eigenvalues and eigenvectors exist
  eigvals <- pca_g$eigenval
  # ensure we have at least 3 eigenvals; pad with NA if missing
  eigen_PC1 <- ifelse(length(eigvals) >= 1, eigvals[1], NA)
  eigen_PC2 <- ifelse(length(eigvals) >= 2, eigvals[2], NA)
  eigen_PC3 <- ifelse(length(eigvals) >= 3, eigvals[3], NA)
  total_var <- sum(eigvals, na.rm = TRUE)
  var_summary <- rbind(var_summary,
                       data.frame(Group = g,
                                  eigen_PC1 = eigen_PC1,
                                  eigen_PC2 = eigen_PC2,
                                  eigen_PC3 = eigen_PC3,
                                  total_variance = total_var,
                                  n_samples = n_samps,
                                  stringsAsFactors = FALSE))
  # store per-sample PC1-3 coordinates for clustering later
  pcs <- pca_g$eigenvect
  # ensure columns exist
  pc1 <- pcs[,1]
  pc2 <- if (ncol(pcs) >= 2) pcs[,2] else rep(0, length(pc1))
  pc3 <- if (ncol(pcs) >= 3) pcs[,3] else rep(0, length(pc1))
  temp <- data.frame(sample = pca_g$sample.id,
                     PC1 = pc1, PC2 = pc2, PC3 = pc3,
                     Group = g, stringsAsFactors = FALSE)
  pca_per_group_samples <- rbind(pca_per_group_samples, temp)

  # optional plot per group (save small figure)
  if (SAVE_PLOTS) {
    if (!dir.exists(OUTDIR)) dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
    pfn <- file.path(OUTDIR, paste0("groupPCA_", g, "_PC1_PC2.png"))
    png(pfn, width = 2000, height = 1600, res = 200)
    plot(temp$PC1, temp$PC2, pch = 19, main = paste0("Group ", g, " PCA (n=", n_samps, ")"),
         xlab = paste0("PC1 (eigen=", signif(eigen_PC1,3), ")"),
         ylab = paste0("PC2 (eigen=", signif(eigen_PC2,3), ")"))
    text(temp$PC1, temp$PC2, labels = temp$sample, cex = 0.6, pos = 3)
    dev.off()
  }
}

OUTDIR <- sub("/+$", "", OUTDIR)
if (!dir.exists(OUTDIR)) dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

subgroup_var_summary <- data.frame(
  Group = character(),
  k = integer(),
  cluster = integer(),
  n_samples = integer(),
  eigen_PC1 = numeric(),
  eigen_PC2 = numeric(),
  eigen_PC3 = numeric(),
  total_variance = numeric(),
  group_total_variance = numeric(),
  stringsAsFactors = FALSE
)

var_summary2 <- subset(var_summary, total_variance > 10 & n_samples > 2)
groups_to_split <- var_summary2$Group
cat("Splitting ALL groups:", paste(groups_to_split, collapse = ", "), "\n")

# collect all plot_df across groups to write a single CSV later
all_plot_dfs <- list()

for (g in groups_to_split) {
  cat("\nProcessing Group:", g, "\n")

  dfg_base <- pca_per_group_samples %>% filter(Group == g)
  if (nrow(dfg_base) < 2) {
    message(sprintf("⚠️ Skipping Group %s: only %d sample(s)", g, nrow(dfg_base)))
    next
  }

  pcs_mat <- as.matrix(dfg_base %>% select(PC1, PC2, PC3))
  stopifnot(all(is.finite(pcs_mat)))
  group_total_var <- var_summary %>% filter(Group == g) %>% pull(total_variance)

  # collect per-k assignments for THIS group only
  plot_df <- dfg_base[0, c("sample","Group","PC1","PC2","PC3")]
  plot_df$k <- integer(0)
  plot_df$cluster <- integer(0)

  for (k in 2:6) {
    if (nrow(dfg_base) < k) {
      message(sprintf("⚠️ Skipping k=%d for Group %s: n=%d too small", k, g, nrow(dfg_base)))
      next
    }

    cat("  k =", k, " -> running kmeans\n")
    set.seed(1)
    km <- kmeans(pcs_mat, centers = k, nstart = 50)

    dfg_k <- dfg_base %>%
      mutate(k = k, cluster = km$cluster) %>%
      select(sample, Group, k, cluster, PC1, PC2, PC3)

    # append to this group's plotting df
    plot_df <- bind_rows(plot_df, dfg_k)

    # PCA-on-genotypes summaries per cluster
    for (cl in 1:k) {
      samples.cluster <- dfg_k$sample[dfg_k$cluster == cl]
      n_cl <- length(samples.cluster)
      cat("    cluster", cl, "n=", n_cl, "\n")

      if (n_cl >= 2) {
        pca_sub <- tryCatch(
          snpgdsPCA(gds, sample.id = samples.cluster, snp.id = snp.id,
                    num.thread = 1, autosome.only = FALSE),
          error = function(e) NULL
        )
        if (!is.null(pca_sub)) {
          eigvals <- pca_sub$eigenval
          eigen1 <- ifelse(length(eigvals) >= 1, eigvals[1], NA)
          eigen2 <- ifelse(length(eigvals) >= 2, eigvals[2], NA)
          eigen3 <- ifelse(length(eigvals) >= 3, eigvals[3], NA)
          total_var <- sum(eigvals, na.rm = TRUE)
        } else {
          eigen1 <- NA; eigen2 <- NA; eigen3 <- NA; total_var <- n_cl
        }
      } else {
        eigen1 <- NA; eigen2 <- NA; eigen3 <- NA; total_var <- n_cl
      }

      subgroup_var_summary <- rbind(subgroup_var_summary, data.frame(
        Group = g, k = k, cluster = cl, n_samples = n_cl,
        eigen_PC1 = eigen1, eigen_PC2 = eigen2, eigen_PC3 = eigen3,
        total_variance = total_var, group_total_variance = group_total_var,
        stringsAsFactors = FALSE
      ))
    }
  } # end k loop

  # Save a per-group faceted plot (optional) reusing plot_df just computed
  if (SAVE_PLOTS && nrow(plot_df) > 0) {
    plot_df$k <- factor(plot_df$k, levels = sort(unique(plot_df$k)))
    p <- ggplot(plot_df, aes(PC1, PC2, color = as.factor(cluster))) +
      geom_point(size = 2, alpha = 0.9) +
      facet_wrap(~ k, ncol = 3, scales = "fixed") +
      coord_equal() +
      labs(
        title = paste0("Group ", g, " - PC1 vs PC2 with k-means clusters"),
        subtitle = "Clustering on Group PC1–PC3; points are Group PCA PC1–PC2",
        x = "PC1 (Group PCA)", y = "PC2 (Group PCA)", color = "Cluster"
      ) +
      theme_bw(base_size = 12) +
      theme(legend.position = "bottom",
            panel.grid.minor = element_blank(),
            plot.title = element_text(face = "bold"))
    pfn <- file.path(OUTDIR, paste0("Group_", g, "_PC1PC2_k2to6.png"))
    message("Saving -> ", pfn); ggsave(pfn, plot = p, width = 12, height = 8, dpi = 300)
    message("Saved: ", pfn)
  }

  # ---- COLLECT THIS GROUP'S plot_df FOR THE MASTER CSV ----
  if (nrow(plot_df) > 0) {
    # ensure stable column order and types
    plot_df <- plot_df %>%
      mutate(Group = as.character(Group),
             k = as.integer(k),
             cluster = as.integer(cluster)) %>%
      select(Group, sample, k, cluster, PC1, PC2, PC3)
    all_plot_dfs[[as.character(g)]] <- plot_df
  } else {
    message(sprintf("[Group %s] no clustering rows to collect", g))
  }
}

# ---- WRITE THE SINGLE BOUND CSV ----
if (length(all_plot_dfs) > 0) {
  all_plot_df <- bind_rows(all_plot_dfs)
  out_csv <- file.path(OUTDIR, "Group_kmeans_assignments_PC1PC2PC3.csv")
  write.csv(all_plot_df, out_csv, row.names = FALSE)
  cat("Wrote combined cluster assignments CSV ->", out_csv, "\n")
} else {
  cat("No plot_df content collected. CSV not written.\n")
}

# Optional: also write subgroup variance summary
subgroup_csv <- file.path(OUTDIR, "Subgroup_PCA_variance_summary.csv")
write.csv(subgroup_var_summary, subgroup_csv, row.names = FALSE)
cat("Wrote subgroup PCA variance summary ->", subgroup_csv, "\n")

write.csv(var_summary, file = file.path(OUTDIR, "group_eigen_summary.csv"), row.names = FALSE)
write.csv(subgroup_var_summary, file = file.path(OUTDIR, "subgroup_eigen_summary.csv"), row.names = FALSE)



