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
het_data <- read.table("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/Gilmer_het_10kb.txt",
                       header = FALSE, stringsAsFactors = FALSE)

colnames(het_data) <- c("contig", "win_start", "win_end", "sample", "het_prop")

# List all files











het_merged <- het_data %>%
  left_join(genomic_types, by = c("sample" = "CloneA")) %>%
  left_join(depths, by = c("sample" = "sampleId"))






# Plot
pdf("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/het_window_Gilmer_new_10000.pdf",
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
      title = paste0("Windowed heterozygosity (10kb): ", s,
                     " | Group: ", group, " | meanDepth: ", depth),
      x = "Genomic position (window start)",
      y = "Observed heterozygosity (Ho)"
    ) +
    ylim(0, 0.05) +
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
      title = "Windowed heterozygosity (10kb)",
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
pdf("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/het_window_Gilmer_10000_overlap.pdf",
    width = 100, height = 10)
p2
dev.off()




  p3 <- ggplot(subset(het_merged, meanDepth > 5), aes(x = win_start, y = het_prop, group= sample, col=meanDepth)) +
    geom_line() +
    facet_grid(~contig, scales = "free_x") +
    labs(
      title = "Windowed heterozygosity (10kb)",
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
pdf("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/het_window_Gilmer_10000_overlap_depth5.pdf",
    width = 100, height = 10)
p3
dev.off()













# Plot
pdf("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/het_window_Gilmer_new_10000_depth5.pdf",
    width = 25, height = 2)

for(s in unique(subset(het_merged, meanDepth > 5)$sample)){
  sample_info <- het_merged %>% filter(sample == s) %>% slice(1)
  group <- sample_info$Group
  depth <- round(sample_info$meanDepth, 2)

  p <- ggplot(subset(het_merged, sample == s),
              aes(x = win_start, y = het_prop)) +
    geom_point(size = 0.4, alpha = 0.6, color = "darkgreen") +
    facet_grid(~contig, scales = "free_x") +
    labs(
      title = paste0("Windowed heterozygosity (10kb): ", s,
                     " | Group: ", group, " | meanDepth: ", depth),
      x = "Genomic position (window start)",
      y = "Observed heterozygosity (Ho)"
    ) +
    ylim(0, 0.05) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.minor = element_blank()
    )
  print(p)
}

dev.off()













library(dplyr)
library(tidyr)

# example dataframe = df

# -----------------------------
# 1. Wide format for meanDepth
# -----------------------------
wide_depth <- het_merged %>%
  select(sample, contig, win_start, meanDepth) %>%
  unite("window", contig, win_start, sep="_") %>%   # make unique window ID
  pivot_wider(names_from = window, values_from = meanDepth)

# Now `wide_depth` has one row per sample, columns = each window meanDepth

# -----------------------------
# 2. Wide format for het_prop
# -----------------------------
wide_het <- het_merged %>%
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

# ----------------------------
# 2. Merge with full metadata
# ----------------------------
# make sure both have same sample IDs
merged <- left_join(pca_df, depths, by = "sampleId")

# ----------------------------
# 3. Variance explained labels
# ----------------------------
var_explained <- summary(pca_res)$importance[2,] * 100
pc1_lab <- paste0("PC1 (", round(var_explained[1],1), "%)")
pc2_lab <- paste0("PC2 (", round(var_explained[2],1), "%)")
pc3_lab <- paste0("PC3 (", round(var_explained[3],1), "%)")
pc4_lab <- paste0("PC4 (", round(var_explained[4],1), "%)")

# ----------------------------
# 4. Plot PCA (color by meanDepth)
# ----------------------------


png("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/pca_het.png",
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





png("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/pca_het_2_3.png",
    width = 500, height = 500)

ggplot(merged, aes(x = PC2, y = PC3, color = meanDepth, label = sampleId)) +
  geom_point(size = 3) +
  scale_color_viridis_c(option = "plasma") +
  theme_minimal() +
  labs(title = "PCA of het_prop per window",
       x = pc2_lab,
       y = pc3_lab,
       color = "Mean Depth")

dev.off()





png("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/pca_het_3_4.png",
    width = 500, height = 500)

ggplot(merged, aes(x = PC3, y = PC4, color = meanDepth, label = sampleId)) +
  geom_point(size = 3) +
  scale_color_viridis_c(option = "plasma") +
  theme_minimal() +
  labs(title = "PCA of het_prop per window",
       x = pc3_lab,
       y = pc4_lab,
       color = "Mean Depth")

dev.off()