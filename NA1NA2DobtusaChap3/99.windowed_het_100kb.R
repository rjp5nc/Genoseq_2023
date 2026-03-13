#ijob -A berglandlab -c10 -p standard --mem=40G
#module load gcc/11.4.0  openmpi/4.1.4 icu R/4.3.1
#R


library(tidyverse)
library(gridExtra)
library(ggplot2)
# Read the vcftools --het output
genomic_types <- read.csv("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/genomic_types.csv")
depths <-   read.csv("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/sampleStats_US_obtusa.csv")

metadata_with_clone <- read.csv("/project/berglandlab/Robert/UKSequencing2022_2024/old_stuff/metadata_with_clone.csv", header = TRUE)
metadata_with_clone$clone <- trimws(metadata_with_clone$clone)
metadata_with_clone <- subset(metadata_with_clone, clone !="Blank")
metadata_with_clone <- subset(metadata_with_clone, clone !="BLANK")


clones_to_update <- c(
  "Rockpool1_D12", "Rockpool1_D8",
  "Rockpool1_C1", "Rockpool1_C2", "Rockpool1_C3",
  "Rockpool2_A1", "Rockpool2_A10", "Rockpool2_A11", "Rockpool2_A12",
  "Rockpool2_A2", "Rockpool2_A3", "Rockpool2_A4", "Rockpool2_A5",
  "Rockpool2_A6", "Rockpool2_A7", "Rockpool2_A8", "Rockpool2_A9",
  "Rockpool2_B1", "Rockpool2_B10", "Rockpool2_B11", "Rockpool2_B12",
  "Rockpool2_B2", "Rockpool2_B3", "Rockpool2_B4", "Rockpool2_B5",
  "Rockpool2_B6", "Rockpool2_B7", "Rockpool2_B8", "Rockpool2_B9",
  "Rockpool2_C2"
)

# add "_2" to Group where CloneA is in the list
genomic_types$Group[genomic_types$CloneA %in% clones_to_update] <-
  paste0(genomic_types$Group[genomic_types$CloneA %in% clones_to_update], "_2")


write.csv(genomic_types, "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/genomic_types_v2.csv")





# Folder where all per-contig het files are
het_data <- read.table("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/only2_het_100kb.txt",
                       header = FALSE, stringsAsFactors = FALSE)

colnames(het_data) <- c("contig", "win_start", "win_end", "sample", "het_prop")

# List all files


het_merged <- het_data %>%
  left_join(genomic_types, by = c("sample" = "CloneA")) %>%
  left_join(depths, by = c("sample" = "sampleId")) %>%
  left_join(metadata_with_clone, by = c("sample" = "Well"))





# Plot
pdf("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/het_window_100000.pdf",
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
      title = paste0("Windowed heterozygosity (100kb): ", s,
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



het_data <- read.table("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/only2_het_10kb_all.txt",
                       header = FALSE, stringsAsFactors = FALSE)

colnames(het_data) <- c("contig", "win_start", "win_end", "sample", "het_prop")

# List all files


het_merged <- het_data %>%
  left_join(genomic_types, by = c("sample" = "CloneA")) %>%
  left_join(depths, by = c("sample" = "sampleId")) %>%
  left_join(metadata_with_clone, by = c("sample" = "Well"))


# Plot
pdf("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/het_window_10000_AC.pdf",
    width = 25, height = 2)

for(s in unique(subset(het_merged, Group == "AC")$sample)){
  sample_info <- het_merged %>% filter(sample == s) %>% slice(1)
  group <- sample_info$Group
  depth <- round(sample_info$meanDepth, 2)

  p <- ggplot(subset(het_merged, sample == s),
              aes(x = win_start, y = het_prop/sqrt(depth))) +
    geom_point(size = 0.4, alpha = 0.6, color = "darkgreen") +
    facet_grid(~contig, scales = "free_x") +
    labs(
      title = paste0("Windowed heterozygosity (10kb): ", s,
                     " | Group: ", group, " | meanDepth: ", depth),
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
  print(p)
}

dev.off()





# Plot
pdf("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/het_window_10000_G_avg.pdf",
    width = 25, height = 2)

for(s in unique(subset(het_merged, Group == "G")$sample)){
  sample_info <- het_merged %>% filter(sample == s) %>% slice(1)
  group <- sample_info$Group
  depth <- round(sample_info$meanDepth, 2)

  p <- ggplot(subset(het_merged, sample == s),
              aes(x = win_start, y = het_prop/sqrt(depth))) +
    geom_point(size = 0.4, alpha = 0.6, color = "darkgreen") +
    facet_grid(~contig, scales = "free_x") +
    labs(
      title = paste0("Windowed heterozygosity (10kb): ", s,
                     " | Group: ", group, " | meanDepth: ", depth),
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
  print(p)
}

dev.off()






# Plot
pdf("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/het_window_10000_G.pdf",
    width = 25, height = 2)

for(s in unique(subset(het_merged, Group == "G")$sample)){
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



# Plot
pdf("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/het_window_10000_M.pdf",
    width = 25, height = 2)

for(s in unique(subset(het_merged, Group == "M")$sample)){
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











het_data <- read.table("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/only2_het_100kb_all.txt",
                       header = FALSE, stringsAsFactors = FALSE)

colnames(het_data) <- c("contig", "win_start", "win_end", "sample", "het_prop")

# List all files


het_merged <- het_data %>%
  left_join(genomic_types, by = c("sample" = "CloneA")) %>%
  left_join(depths, by = c("sample" = "sampleId")) %>%
  left_join(metadata_with_clone, by = c("sample" = "Well"))


het_mergedgilmer <- subset(het_merged, accuratelocation == "Gilmer")

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

png("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/pca_Gilmer_A_het.png",
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






het_mergedG <- subset(het_merged, Group == "G")

wide_het <- het_mergedG %>%
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

png("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/pca_G_het.png",
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

library(data.table)

some <- setDT(pca_df, keep.rownames = "Sample")

merged <- left_join(pca_df, depths, by = "sampleId")



merged2 <- pca_df %>%
  left_join(depths, by = c("sampleId" = "sampleId")) %>%
  left_join(metadata_with_clone, by = c("sample" = "Well"))



var_explained <- summary(pca_res)$importance[2,] * 100
pc1_lab <- paste0("PC1 (", round(var_explained[1],1), "%)")
pc2_lab <- paste0("PC2 (", round(var_explained[2],1), "%)")
pc3_lab <- paste0("PC3 (", round(var_explained[3],1), "%)")
pc4_lab <- paste0("PC4 (", round(var_explained[4],1), "%)")

png("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/pca_all_het.png",
    width = 5000, height = 5000)

ggplot(merged, aes(x = PC1, y = PC2, color = meanDepth, label = sampleId)) +
  geom_point(size = 3) +
  scale_color_viridis_c(option = "plasma") +
  theme_minimal() +
  facet_wrap("Group") +
  labs(title = "PCA of het_prop per window",
       x = pc1_lab,
       y = pc2_lab,
       color = "Mean Depth")

dev.off()










# Plot
pdf("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/het_window_100k_A_all.pdf",
    width = 50, height = 10)
  p <- ggplot(subset(het_merged, Group == "A" & meanDepth >= 5),
              aes(x = win_start, y = het_prop/sqrt(meanDepth), group =sample, col = meanDepth)) +
    geom_line() +
    facet_grid(date~contig, scales = "free_x") +
    labs(
      title = paste0("Windowed heterozygosity (100kb): "),
      x = "Genomic position (window start)",
      y = "Observed heterozygosity (Ho)"
    ) +
#    ylim(0, 0.005) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.minor = element_blank()
    )
  print(p)
dev.off()





# Plot
pdf("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/het_window_100k_F.pdf",
    width = 25, height = 2)

for(s in unique(subset(het_merged, Group == "F")$sample)){
  sample_info <- het_merged %>% filter(sample == s) %>% slice(1)
  group <- sample_info$Group
  depth <- round(sample_info$meanDepth, 2)

  p <- ggplot(subset(het_merged, sample == s),
              aes(x = win_start, y = het_prop)) +
    geom_point(size = 0.4, alpha = 0.6, color = "darkgreen") +
    facet_grid(~contig, scales = "free_x") +
    labs(
      title = paste0("Windowed heterozygosity (100kb): ", s,
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




# Plot
pdf("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/het_window_100k_T.pdf",
    width = 25, height = 2)

for(s in unique(subset(het_merged, Group == "T")$sample)){
  sample_info <- het_merged %>% filter(sample == s) %>% slice(1)
  group <- sample_info$Group
  depth <- round(sample_info$meanDepth, 2)

  p <- ggplot(subset(het_merged, sample == s),
              aes(x = win_start, y = het_prop)) +
    geom_point(size = 0.4, alpha = 0.6, color = "darkgreen") +
    facet_grid(~contig, scales = "free_x") +
    labs(
      title = paste0("Windowed heterozygosity (100kb): ", s,
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

#unsupervised clustering 
#use as many dimensions of the pca as you want/can
#






for (g in unique(het_merged$Group)) {
  # Subset data for this group
  group_data <- subset(het_merged, Group == g)

  # Open PDF for this group
  pdf(paste0("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/het_window_100k_", g, ".pdf"),
      width = 25, height = 2)

  # Loop through each sample in this group
  for (s in unique(group_data$sample)) {
    sample_info <- group_data %>% filter(sample == s) %>% slice(1)
    depth <- round(sample_info$meanDepth, 2)

    p <- ggplot(filter(group_data, sample == s),
                aes(x = win_start, y = het_prop)) +
      geom_point(size = 0.4, alpha = 0.6, color = "darkgreen") +
      facet_grid(~contig, scales = "free_x") +
      labs(
        title = paste0("Windowed heterozygosity (100kb): ", s,
                       " | Group: ", g, " | meanDepth: ", depth),
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
}





























# Plot
pdf("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/het_window_100k_G.pdf",
    width = 25, height = 2)

for(s in unique(subset(het_merged, Group == "G")$sample)){
  sample_info <- het_merged %>% filter(sample == s) %>% slice(1)
  group <- sample_info$Group
  depth <- round(sample_info$meanDepth, 2)

  p <- ggplot(subset(het_merged, sample == s),
              aes(x = win_start, y = het_prop)) +
    geom_point(size = 0.4, alpha = 0.6, color = "darkgreen") +
    facet_grid(~contig, scales = "free_x") +
    labs(
      title = paste0("Windowed heterozygosity (100kb): ", s,
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


# Plot
pdf("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/het_window_100k_AB_E.pdf",
    width = 25, height = 2)

for(s in unique(subset(het_merged, Group == "AB" |  Group == "E")$sample)){
  sample_info <- het_merged %>% filter(sample == s) %>% slice(1)
  group <- sample_info$Group
  depth <- round(sample_info$meanDepth, 2)

  p <- ggplot(subset(het_merged, sample == s),
              aes(x = win_start, y = het_prop)) +
    geom_point(size = 0.4, alpha = 0.6, color = "darkgreen") +
    facet_grid(~contig, scales = "free_x") +
    labs(
      title = paste0("Windowed heterozygosity (100kb): ", s,
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



# Plot
pdf("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/het_window_100k_AC_E_F.pdf",
    width = 25, height = 2)

for(s in unique(subset(het_merged, Group == "AC" |  Group == "E"|  Group == "F")$sample)){
  sample_info <- het_merged %>% filter(sample == s) %>% slice(1)
  group <- sample_info$Group
  depth <- round(sample_info$meanDepth, 2)

  p <- ggplot(subset(het_merged, sample == s),
              aes(x = win_start, y = het_prop)) +
    geom_point(size = 0.4, alpha = 0.6, color = "darkgreen") +
    facet_grid(~contig, scales = "free_x") +
    labs(
      title = paste0("Windowed heterozygosity (100kb): ", s,
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





# Plot
pdf("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/het_window_100k_G_all.pdf",
    width = 50, height = 10)
  p <- ggplot(subset(het_merged, Group == "G"),
              aes(x = win_start, y = het_prop/sqrt(depth), group =sample)) +
    geom_line() +
    facet_grid(~contig, scales = "free_x") +
    labs(
      title = paste0("Windowed heterozygosity (100kb): "),
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
dev.off()












het_data <- read.table("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/only2_het_1000kb_all.txt",
                       header = FALSE, stringsAsFactors = FALSE)

colnames(het_data) <- c("contig", "win_start", "win_end", "sample", "het_prop")

# List all files


het_merged <- het_data %>%
  left_join(genomic_types, by = c("sample" = "CloneA")) %>%
  left_join(depths, by = c("sample" = "sampleId")) %>%
  left_join(metadata_with_clone, by = c("sample" = "Well"))


  





# Plot
pdf("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/het_window_1000k_G.pdf",
    width = 25, height = 2)

for(s in unique(subset(het_merged, Group == "G")$sample)){
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


# Plot
pdf("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/het_window_1000k_G_all.pdf",
    width = 50, height = 10)
  p <- ggplot(subset(het_merged, Group == "G" & meanDepth >= 5),
              aes(x = win_start, y = het_prop/sqrt(meanDepth), group =sample, col = meanDepth)) +
    geom_line() +
    facet_grid(~contig, scales = "free_x") +
    labs(
      title = paste0("Windowed heterozygosity (1000kb): "),
      x = "Genomic position (window start)",
      y = "Observed heterozygosity (Ho)"
    ) +
#    ylim(0, 0.005) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.minor = element_blank()
    )
  print(p)
dev.off()








# Plot
pdf("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/het_window_1000k_A_all.pdf",
    width = 50, height = 10)
  p <- ggplot(subset(het_merged, Group == "A" & meanDepth >= 5),
              aes(x = win_start, y = het_prop/sqrt(meanDepth), group =sample, col = meanDepth)) +
    geom_line() +
    facet_grid(date~contig, scales = "free_x") +
    labs(
      title = paste0("Windowed heterozygosity (1000kb): "),
      x = "Genomic position (window start)",
      y = "Observed heterozygosity (Ho)"
    ) +
#    ylim(0, 0.005) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.minor = element_blank()
    )
  print(p)
dev.off()