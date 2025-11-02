#ijob -A berglandlab -c2 -p standard --mem=40G
#module load gcc/11.4.0  openmpi/4.1.4 icu R/4.3.1
#R


library(tidyverse)
library(gridExtra)
library(ggplot2)
# Read the vcftools --het output
genomic_types <- read.csv("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/genomic_types.csv")
depths <-   read.csv("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/sampleStats_US_obtusa.csv")
# Folder where all per-contig het files are
het_data <- read.table("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/Gilmer_het_1000kb.txt",
                       header = FALSE, stringsAsFactors = FALSE)

colnames(het_data) <- c("contig", "win_start", "win_end", "sample", "het_prop")

# List all files











het_merged <- het_data %>%
  left_join(genomic_types, by = c("sample" = "CloneA")) %>%
  left_join(depths, by = c("sample" = "sampleId"))






# Plot
pdf("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/het_window_Gilmer_new_1000kb.pdf",
    width = 25, height = 2)

for(s in unique(het_merged$sample)){
  sample_info <- het_merged %>% filter(sample == s) %>% slice(1)
  group <- sample_info$Group
  depth <- round(sample_info$meanDepth, 2)

  p <- ggplot(subset(het_merged, sample == s),
              aes(x = win_start, y = het_prop)) +
    geom_point(size = 0.4, alpha = 0.6, color = "darkgreen") +
    facet_grid(~contig, scales = "free_x") +
    labs(
      title = paste0("Windowed heterozygosity (1000kb): ", s,
                     " | Group: ", group, " | meanDepth: ", depth),
      x = "Genomic position (window start)",
      y = "Observed heterozygosity (Ho)"
    ) +
    ylim(0, 0.01) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.minor = element_blank()
    )
  print(p)
}

dev.off()







  p2 <- ggplot(het_merged, aes(x = win_start, y = het_prop, group= sample, col=meanDepth)) +
    geom_line() +
    facet_grid(~contig, scales = "free_x") +
    labs(
      title = "Windowed heterozygosity (1000kb)",
      x = "Genomic position (window start)",
      y = "Observed heterozygosity (Ho)"
    ) +
    ylim(0, 0.02) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.minor = element_blank()
    )



# Plot
pdf("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/het_window_Gilmer_1000000_overlap.pdf",
    width = 100, height = 10)
p2
dev.off()




  p3 <- ggplot(subset(het_merged, meanDepth > 3), aes(x = win_start, y = het_prop, group= sample, col=meanDepth)) +
    geom_line() +
    facet_grid(~contig, scales = "free_x") +
    labs(
      title = "Windowed heterozygosity (1000kb)",
      x = "Genomic position (window start)",
      y = "Observed heterozygosity (Ho)"
    ) +
    ylim(0, 0.02) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.minor = element_blank()
    )



# Plot
pdf("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/het_window_Gilmer_1000000_overlap_depth3.pdf",
    width = 100, height = 10)
p3
dev.off()















het_data <- read.table("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/Gilmer_het_1000kb.txt",
                       header = FALSE, stringsAsFactors = FALSE)

colnames(het_data) <- c("contig", "win_start", "win_end", "sample", "het_prop")

# List all files


het_merged <- het_data %>%
  left_join(genomic_types, by = c("sample" = "CloneA")) %>%
  left_join(depths, by = c("sample" = "sampleId")) 


het_mergedgilmer <- het_merged

wide_het <- het_mergedgilmer %>%
  select(sample, contig, win_start, het_prop) %>%
  unite("window", contig, win_start, sep="_") %>%
  pivot_wider(names_from = window, values_from = het_prop)

# -----------------------------
# 3. PCA on het_prop
# -----------------------------
# Drop "sample" column to run PCA
mat <- wide_het %>%
  column_to_rownames("sample")

# replace NAs with 0 or mean (depending on preference)
mat[is.na(mat)] <- 0

# Remove constant columns (sd = 0)
mat_clean <- mat[, apply(mat, 2, sd) != 0]

# Run PCA
pca_res <- prcomp(mat_clean, scale. = TRUE)

# PCA results
summary(pca_res)

pca_df <- as.data.frame(pca_res$x)
pca_df$sampleId <- rownames(pca_df)

merged <- left_join(pca_df, depths, by = "sampleId")

var_explained <- summary(pca_res)$importance[2,] * 100
pc1_lab <- paste0("PC1 (", round(var_explained[1],1), "%)")
pc2_lab <- paste0("PC2 (", round(var_explained[2],1), "%)")
pc3_lab <- paste0("PC3 (", round(var_explained[3],1), "%)")
pc4_lab <- paste0("PC4 (", round(var_explained[4],1), "%)")

png("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/pca_Gilmer_A_1000kb_het.png",
    width = 500, height = 500)

ggplot(merged, aes(x = PC1, y = PC2, color = meanDepth, label = sampleId)) +
  geom_point(size = 3) +
  scale_color_viridis_c(option = "plasma") +
  theme_minimal() +
  labs(title = "PCA of het_prop per window",
       x = pc1_lab,
       y = pc2_lab,
       color = "Mean Depth")

dev.off()



pc1filt<-subset(merged, PC2 <= -10)$sampleId




















ids_sub <- merged %>%
  dplyr::filter(PC2 >= -10) %>%
  dplyr::pull(sampleId) %>%
  unique()

# --- subset your original feature matrix used for PCA ---
# Expectation: rownames(het_mat) are sample IDs used in the original PCA
sub_mat <- mat_clean[rownames(mat_clean) %in% ids_sub, , drop = FALSE]

# Safety: drop columns with zero variance and optionally impute NAs
nzv <- apply(sub_mat, 2, function(x) stats::var(x, na.rm = TRUE) > 0)
sub_mat <- sub_mat[, nzv, drop = FALSE]

# Simple mean impute per column if needed
na_cols <- which(colSums(is.na(sub_mat)) > 0)
if (length(na_cols) > 0) {
  for (j in na_cols) {
    m <- mean(sub_mat[, j], na.rm = TRUE)
    sub_mat[is.na(sub_mat[, j]), j] <- m
  }
}

# --- run PCA on the subset ---
pca_res_sub <- prcomp(sub_mat, center = TRUE, scale. = TRUE)

# Scores dataframe
pca_df_sub <- as.data.frame(pca_res_sub$x)
pca_df_sub$sampleId <- rownames(pca_res_sub$x)

# Join extra info for coloring
merged_sub <- dplyr::left_join(pca_df_sub, depths, by = "sampleId")

# Variance explained labels
var_explained_sub <- (pca_res_sub$sdev^2) / sum(pca_res_sub$sdev^2) * 100
pc1_lab_sub <- paste0("PC1 (", round(var_explained_sub[1], 1), "%)")
pc2_lab_sub <- paste0("PC2 (", round(var_explained_sub[2], 1), "%)")

# --- plot ---
png("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/pca_Gilmer_A_1000kb_het_subset_PC2le-10.png",
    width = 500, height = 500)

ggplot(merged_sub, aes(x = PC1, y = PC2, color = meanDepth, label = sampleId)) +
  geom_point(size = 3) +
  scale_color_viridis_c(option = "plasma") +
  theme_minimal() +
  labs(title = "PCA of het_prop per window (PC2 ≤ -10 subset)",
       x = pc1_lab_sub,
       y = pc2_lab_sub,
       color = "Mean Depth")

dev.off()


merged_sub


subset(merged, PC2 <= -10)$sampleId


pc2filt<-subset(merged_sub, PC2 <= -1)$sampleId



ids_sub <- merged_sub %>%
  dplyr::filter(PC2 >= -1) %>%
  dplyr::pull(sampleId) %>%
  unique()

# --- subset your original feature matrix used for PCA ---
# Expectation: rownames(het_mat) are sample IDs used in the original PCA
sub_mat <- mat_clean[rownames(mat_clean) %in% ids_sub, , drop = FALSE]

# Safety: drop columns with zero variance and optionally impute NAs
nzv <- apply(sub_mat, 2, function(x) stats::var(x, na.rm = TRUE) > 0)
sub_mat <- sub_mat[, nzv, drop = FALSE]

# Simple mean impute per column if needed
na_cols <- which(colSums(is.na(sub_mat)) > 0)
if (length(na_cols) > 0) {
  for (j in na_cols) {
    m <- mean(sub_mat[, j], na.rm = TRUE)
    sub_mat[is.na(sub_mat[, j]), j] <- m
  }
}

# --- run PCA on the subset ---
pca_res_sub <- prcomp(sub_mat, center = TRUE, scale. = TRUE)

# Scores dataframe
pca_df_sub <- as.data.frame(pca_res_sub$x)
pca_df_sub$sampleId <- rownames(pca_res_sub$x)

# Join extra info for coloring
merged_sub <- dplyr::left_join(pca_df_sub, depths, by = "sampleId")

# Variance explained labels
var_explained_sub <- (pca_res_sub$sdev^2) / sum(pca_res_sub$sdev^2) * 100
pc1_lab_sub <- paste0("PC1 (", round(var_explained_sub[1], 1), "%)")
pc2_lab_sub <- paste0("PC2 (", round(var_explained_sub[2], 1), "%)")

# --- plot ---
png("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/pca_Gilmer_A_1000kb_het_subset_3.png",
    width = 500, height = 500)

ggplot(merged_sub, aes(x = PC1, y = PC2, color = meanDepth, label = sampleId)) +
  geom_point(size = 3) +
  scale_color_viridis_c(option = "plasma") +
  theme_minimal() +
  labs(title = "PCA of het_prop per window (PC2 ≤ -10 subset)",
       x = pc1_lab_sub,
       y = pc2_lab_sub,
       color = "Mean Depth")

dev.off()

pc3filt<-subset(merged_sub, PC2 <= -1)$sampleId

pc1filt
pc2filt
pc3filt <- merged_sub$sampleId


pc1filt
pc2filt
pc3filt


df_combined <- dplyr::bind_rows(
  data.frame(Well = pc1filt, sub_super = "A_1"),
  data.frame(Well = pc2filt, sub_super = "A_2"),
  data.frame(Well = pc3filt, sub_super = "A_3")
)

# Optional: inspect
head(df_combined)
nrow(df_combined)
