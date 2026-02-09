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
het_data <- read.table("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/only2_het_100kb_all.txt",
                       header = FALSE, stringsAsFactors = FALSE)

colnames(het_data) <- c("contig", "win_start", "win_end", "sample", "het_prop")

# List all files




het_merged <- het_data %>%
  left_join(genomic_types, by = c("sample" = "CloneA")) %>%
  left_join(depths, by = c("sample" = "sampleId"))






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
merged <- left_join(merged, genomic_types, by = c("sampleId" = "CloneA"))

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


png("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/pca_het_all_grouped.png",
    width = 2000, height = 2000)

ggplot(merged, aes(x = PC1, y = PC2, color = Group, label = sampleId)) +
  geom_point(size = 3) +
  theme_bw() + facet_wrap(~Group)+
  labs(title = "PCA of het_prop per window",
       x = pc1_lab,
       y = pc2_lab,
       color = "Group")

dev.off()





library(tidyverse)
library(foreach)
library(doParallel)
library(gridExtra)              # load it each session


ncores <- parallel::detectCores() - 1
cl <- makeCluster(ncores)
registerDoParallel(cl)

# -----------------------------------------------------
# 2. Get list of unique groups
# -----------------------------------------------------
groups <- unique(het_merged$Group)

library(foreach)
library(ggplot2)
library(dplyr)
library(tidyr)

foreach(g = unique(het_merged$Group)) %do% {

  # Subset for this group
  wide_het <- het_merged %>%
    filter(Group == g) %>%
    select(sample, contig, win_start, het_prop) %>%
    unite("window", contig, win_start, sep="_") %>%
    pivot_wider(names_from = window, values_from = het_prop)

  mat <- wide_het %>% column_to_rownames("sample")
  mat[is.na(mat)] <- 0
  mat_clean <- mat[, apply(mat, 2, sd) != 0]

  pca_res <- prcomp(mat_clean, scale. = TRUE)
  pca_df <- as.data.frame(pca_res$x)
  pca_df$sampleId <- rownames(pca_df)

  merged <- left_join(pca_df, depths, by = "sampleId") %>%
    left_join(., genomic_types, by = c("sampleId" = "CloneA"))

  var_explained <- summary(pca_res)$importance[2,] * 100
  pc1_lab <- paste0("PC1 (", round(var_explained[1],1), "%)")
  pc2_lab <- paste0("PC2 (", round(var_explained[2],1), "%)")
  pc3_lab <- if(length(var_explained) >= 3) paste0("PC3 (", round(var_explained[3],1), "%)") else NULL

  plots <- list()

  # PC1 vs PC2
  p1 <- ggplot(merged, aes(x = PC1, y = PC2, color = Group, label = sampleId)) +
    geom_point(size = 3) +
    theme_bw() +
    labs(title = paste("PCA of het_prop per window (Group:", g, ") — PC1 vs PC2"),
         x = pc1_lab, y = pc2_lab, color = "Group")
  plots[[1]] <- p1

  # PC2 vs PC3 (only if exists)
  if ("PC3" %in% colnames(merged)) {
    p2 <- ggplot(merged, aes(x = PC2, y = PC3, color = Group, label = sampleId)) +
      geom_point(size = 3) +
      theme_bw() +
      labs(title = paste("PCA of het_prop per window (Group:", g, ") — PC2 vs PC3"),
           x = pc2_lab, y = pc3_lab, color = "Group")
    plots[[2]] <- p2
  }

  # Save plot safely
  out_path <- paste0("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/pca_het_", g, ".png")
  png(out_path, width = 2000, height = 2000, res = 300)

  if (length(plots) == 2) {
    if (requireNamespace("patchwork", quietly = TRUE)) {
      library(patchwork)
      print(plots[[1]] / plots[[2]])  # stack vertically
    } else if (requireNamespace("gridExtra", quietly = TRUE)) {
      library(gridExtra)
      do.call(grid.arrange, c(plots, ncol = 1))
    } else {
      print(plots[[1]])
      print(plots[[2]])
    }
  } else {
    print(plots[[1]])
  }

  dev.off()

  message(paste("Finished PCA for group:", g))
}

# -----------------------------------------------------
# 4. Stop cluster
# -----------------------------------------------------
stopCluster(cl)































































png("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/pca_het_all.png",
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



png("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/pca_het_all_subsets.png",
    width = 3000, height = 3000, res = 600)

ggplot(subset(merged, Group %in% c("C", "G", "H", "N")), aes(x = PC1, y = PC2, color = Group, label = sampleId)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PCA of het_prop per window",
       x = pc1_lab,
       y = pc2_lab,
       color = "Group")

dev.off()





png("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/pca_het_all_grouped_2_3.png",
    width = 2000, height = 2000)

ggplot(merged, aes(x = PC2, y = PC3, color = Group, label = sampleId)) +
  geom_point(size = 3) +
  theme_bw() + facet_wrap(~Group)+
  labs(title = "PCA of het_prop per window",
       x = pc2_lab,
       y = pc3_lab,
       color = "Mean Depth")

dev.off()









k <- 2  # choose number of clusters

mergedF <- subset(merged, Group == "F")



set.seed(123)
library(cluster)

sil_width <- sapply(2:10, function(k) {
  km <- kmeans(mergedF[, 2:5], centers = k, nstart = 25)
  ss <- silhouette(km$cluster, dist(mergedF[, 2:5]))
  mean(ss[, 3])  # average silhouette width
})


png("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/pca_het_all_grouped_2_3_cluster_silwidth.png",
    width = 2000, height = 2000)

plot(2:10, sil_width, type = "b", pch = 19, frame = FALSE,
     xlab = "Number of clusters (k)",
     ylab = "Average silhouette width",
     main = "Silhouette Method for Optimal k")

dev.off()

optimal_k <- which.max(sil_width) + 1  # since we started at k=2
optimal_k

set.seed(123)
km <- kmeans(mergedF[, 2:5], centers = optimal_k, nstart = 25)
mergedF$cluster <- as.factor(km$cluster)



png("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/pca_het_all_grouped_2_3_cluster.png",
    width = 2000, height = 2000)

ggplot(mergedF, aes(PC2, PC3, color = cluster)) +
  geom_point(size = 3) +
  theme_bw() +
  scale_color_viridis_d() +
  labs(title = paste("K-means Clustering (k =", optimal_k, ")"))

dev.off()







merged %>% summarise(
  n_PC1_NA = sum(is.na(PC1)),
  n_PC2_NA = sum(is.na(PC2))
)





merged2 <- merged %>%
  filter(!is.na(Group)) %>%           # remove rows with NA Group
  group_by(Group) %>%
  group_modify(~ {
    df <- .x
    n_ind <- nrow(df)
    
    # Only run kmeans if n_ind > 1 and choose up to 5 clusters
    if (n_ind > 1) {
      n_clusters <- min(5, n_ind - 1)  # kmeans requires n_clusters < n_ind
      df$SubGroup <- kmeans(cbind(df$PC1, df$PC2), centers = n_clusters)$cluster
    } else {
      df$SubGroup <- 1
    }
    
    df
  }) %>%
  ungroup()


png("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/pca_het_all_subgroups.png",
    width = 2000, height = 2000, res = 300)

ggplot(merged2, aes(x = PC1, y = PC2, color = Group, shape = as.factor(SubGroup))) +
  geom_point(size = 3, alpha = 0.8) +                                   # points
  stat_ellipse(aes(group = interaction(Group, SubGroup)), linetype = 2) + # ellipses
  theme_bw(base_size = 16) +
  labs(title = "PCA of heterozygosity (subgroups)",
       x = pc1_lab,
       y = pc2_lab,
       color = "Group",
       shape = "SubGroup") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "right")

dev.off()
















library(dplyr)
library(tidyr)
library(ggplot2)

# 1️⃣ Compute the average het_prop per sample
het_binary <- het_merged %>%
  group_by(sample) %>%
  mutate(sample_mean_het = mean(het_prop, na.rm = TRUE)) %>%
  ungroup() %>%
  # Binarize: 1 if het_prop > sample_mean_het, else 0
  mutate(het_bin = ifelse(het_prop > sample_mean_het, 1, 0))

# 2️⃣ Pivot to wide format: samples as rows, windows as columns
wide_bin <- het_binary %>%
  select(sample, contig, win_start, het_bin) %>%
  unite("window", contig, win_start, sep="_") %>%
  pivot_wider(names_from = window, values_from = het_bin, values_fill = 0)  # fill NA with 0

# 3️⃣ Convert to matrix for PCA
mat_bin <- wide_bin %>% column_to_rownames("sample")
mat_bin <- as.matrix(mat_bin)

# Optional: remove non-variable columns
mat_bin <- mat_bin[, apply(mat_bin, 2, sd) != 0]

# 4️⃣ Run PCA
pca_res <- prcomp(mat_bin, scale. = TRUE)
pca_df <- as.data.frame(pca_res$x)
pca_df$sampleId <- rownames(pca_df)



png("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/pca_het_binary_1_2.png",
    width = 5000, height = 5000, res = 600)

# 5️⃣ Visualize PC1 vs PC2
ggplot(pca_df, aes(x = PC1, y = PC2, label = sampleId)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.5, size = 3) +
  theme_bw() +
  labs(title = "PCA on Binary Het Prop (> sample mean)",
       x = paste0("PC1 (", round(summary(pca_res)$importance[2,1]*100,1), "%)"),
       y = paste0("PC2 (", round(summary(pca_res)$importance[2,2]*100,1), "%)"))

dev.off()






library(dplyr)
library(tidyr)
library(ggplot2)

# 1️⃣ Compute the average het_prop per sample and binarize
het_binary <- het_merged %>%
  group_by(sample) %>%
  mutate(sample_mean_het = mean(het_prop, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(het_bin = ifelse(het_prop > sample_mean_het, 1, 0))

het_binary2 <- as.data.frame(het_binary)
het_binary2G <- subset(het_binary2, Group == "G" & contig =="JAACYE010000001.1")


table_by_sample <- het_binary2G %>%
  group_by(sample) %>%
  summarize(
    n_zeros = sum(het_bin == 0, na.rm = TRUE),
    n_ones  = sum(het_bin == 1, na.rm = TRUE),
    total_windows = n()
  ) %>%
  ungroup()

table_by_sample <- as.data.frame(table_by_sample)
table_by_sample



# 2️⃣ Prepare a dataframe to store PCA results
pca_all_groups <- data.frame()

groups <- unique(het_binary$Group)

for(g in groups){

  # Skip NA group
  if(is.na(g)) next

  # Subset group
  het_grp <- het_binary %>% filter(Group == g)

  # Pivot to wide format
  wide_bin <- het_grp %>%
    select(sample, contig, win_start, het_bin) %>%
    unite("window", contig, win_start, sep="_") %>%
    pivot_wider(names_from = window, values_from = het_bin, values_fill = 0)

  # Skip groups with <2 samples
  if(nrow(wide_bin) < 2) next

  # Convert to matrix and remove non-variable columns
  mat_bin <- wide_bin %>% column_to_rownames("sample") %>% as.matrix()
  mat_bin <- mat_bin[, apply(mat_bin, 2, sd) != 0]

  if(ncol(mat_bin) < 1) next

  # PCA
  pca_res <- prcomp(mat_bin, scale. = TRUE)
  pca_df <- as.data.frame(pca_res$x[, 1:2])  # keep PC1 and PC2
  pca_df$sampleId <- rownames(pca_df)
  pca_df$Group <- g

  # Combine
  pca_all_groups <- rbind(pca_all_groups, pca_df)
}

# 3️⃣ Plot all groups using facet_wrap
png("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/pca_het_binary_facet.png",
    width = 10000, height = 5000, res = 600)

ggplot(pca_all_groups, aes(x = PC1, y = PC2, label = sampleId)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.5, size = 2.5) +
  facet_wrap(~Group) +
#  facet_wrap(~Group, scales = "free") +
  theme_bw() +
  labs(title = "PCA on Binary Het Prop (> sample mean) by Group",
       x = "PC1", y = "PC2")

dev.off()



#PC3

pca_all_groups <- data.frame()
groups <- unique(het_binary$Group)

for(g in groups){
  
  if(is.na(g)) next

  # Subset group
  het_grp <- het_binary %>% filter(Group == g)

  # Pivot to wide format
  wide_bin <- het_grp %>%
    select(sample, contig, win_start, het_bin) %>%
    unite("window", contig, win_start, sep="_") %>%
    pivot_wider(names_from = window, values_from = het_bin, values_fill = 0)

  if(nrow(wide_bin) < 2) next

  # Convert to matrix and remove non-variable columns
  mat_bin <- wide_bin %>% column_to_rownames("sample") %>% as.matrix()
  mat_bin <- mat_bin[, apply(mat_bin, 2, sd) != 0]
  if(ncol(mat_bin) < 1) next

  # PCA
  pca_res <- prcomp(mat_bin, scale. = TRUE)

  # Ensure 3 PCs exist
  n_pcs <- ncol(pca_res$x)
  pc1 <- pca_res$x[,1]
  pc2 <- pca_res$x[,2]
  pc3 <- if(n_pcs >= 3) pca_res$x[,3] else rep(0, nrow(pca_res$x))

  pca_df <- data.frame(
    PC1 = pc1,
    PC2 = pc2,
    PC3 = pc3,
    sampleId = rownames(pca_res$x),
    Group = g
  )

  pca_all_groups <- rbind(pca_all_groups, pca_df)
}


png("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/pca_het_binary_facet_2_3.png",
    width = 10000, height = 5000, res = 600)

ggplot(pca_all_groups, aes(x = PC2, y = PC3, label = sampleId)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.5, size = 2.5) +
  facet_wrap(~Group) +
#  facet_wrap(~Group, scales = "free") +
  theme_bw() +
  labs(title = "PCA on Binary Het Prop (> sample mean) by Group",
       x = "PC2", y = "PC3")

dev.off()



png("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/pca_het_binary_facet_2_3_A.png",
    width = 10000, height = 5000, res = 600)

ggplot(subset(pca_all_groups, Group == "A"), aes(x = PC2, y = PC3, label = sampleId)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.5, size = 2.5) +
  facet_wrap(~Group) +
#  facet_wrap(~Group, scales = "free") +
  theme_bw() +
  labs(title = "PCA on Binary Het Prop (> sample mean) by Group",
       x = "PC2", y = "PC3")

dev.off()




png("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/pca_het_binary_facet_2_3_G.png",
    width = 10000, height = 5000, res = 600)

ggplot(subset(pca_all_groups, Group == "G"), aes(x = PC2, y = PC3, label = sampleId)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.5, size = 2.5) +
  facet_wrap(~Group) +
#  facet_wrap(~Group, scales = "free") +
  theme_bw() +
  labs(title = "PCA on Binary Het Prop (> sample mean) by Group",
       x = "PC2", y = "PC3")

dev.off()





png("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/pca_het_binary_facet_1_2_G.png",
    width = 10000, height = 5000, res = 600)

ggplot(subset(pca_all_groups, Group == "G"), aes(x = PC1, y = PC2, label = sampleId)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.5, size = 2.5) +
  facet_wrap(~Group) +
#  facet_wrap(~Group, scales = "free") +
  theme_bw() +
  labs(title = "PCA on Binary Het Prop (> sample mean) by Group",
       x = "PC1", y = "PC2")

dev.off()



png("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/pca_het_binary_facet_1_2_A.png",
    width = 10000, height = 5000, res = 600)

ggplot(subset(pca_all_groups, Group == "A"), aes(x = PC1, y = PC2, label = sampleId)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.5, size = 2.5) +
  facet_wrap(~Group) +
#  facet_wrap(~Group, scales = "free") +
  theme_bw() +
  labs(title = "PCA on Binary Het Prop (> sample mean) by Group",
       x = "PC1", y = "PC2")

dev.off()




png("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/pca_het_binary_facet_1_2.png",
    width = 10000, height = 5000, res = 600)

ggplot(subset(pca_all_groups), aes(x = PC1, y = PC2, label = sampleId)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.5, size = 2.5) +
  facet_wrap(~Group) +
#  facet_wrap(~Group, scales = "free") +
  theme_bw() +
  labs(title = "PCA on Binary Het Prop (> sample mean) by Group",
       x = "PC1", y = "PC2")

dev.off()


































suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

# --------------------------
# 1) Clean metadata and join
# --------------------------
metadata_with_clone <- read.csv(
  "/project/berglandlab/Robert/UKSequencing2022_2024/old_stuff/metadata_with_clone.csv",
  header = TRUE, stringsAsFactors = FALSE
) %>%
  filter(!(clone %in% c("Blank","BLANK")))

# Optional: sanity-check overlap
# sum(het_merged$sample %in% metadata_with_clone$Well)

# --------------------------
# 2) Per-sample genome-wide heterozygosity
#    (mean across windows; tweak filters if needed)
# --------------------------
# Example quality filters (uncomment/tweak if you want):
# het_merged <- het_merged %>%
#   filter(missingRate <= 0.2, meanDepth >= 5)

het_per_sample <- het_merged %>%
  group_by(sample) %>%
  summarise(
    mean_het   = mean(het_prop, na.rm = TRUE),
    median_het = median(het_prop, na.rm = TRUE),
    n_windows  = sum(!is.na(het_prop)),
    mean_depth = mean(meanDepth, na.rm = TRUE),
    .groups = "drop"
  )

# --------------------------
# 3) Merge with metadata
# --------------------------
merged <- het_per_sample %>%
  left_join(
    metadata_with_clone %>%
      select(Well, date, location, accuratelocation, Plate, clone, sex),
    by = c("sample" = "Well")
  ) %>%
  mutate(
    pool = ifelse(!is.na(accuratelocation) & nzchar(accuratelocation),
                  accuratelocation, location),
    timepoint = as.character(date)  # if 'date' is a year; otherwise adapt parse
  )

# Optional: flag samples that failed to merge
# merged %>% filter(is.na(pool)) %>% distinct(sample)

# --------------------------
# 4) Average per pool per timepoint
# --------------------------
pool_time_summary <- pool_time_summary %>%
  mutate(
    # Try to parse flexible date formats
    parsed_time = suppressWarnings(mdy(timepoint)),           # handles "5/4/2024"
    parsed_time = if_else(is.na(parsed_time),
                          suppressWarnings(ymd(timepoint)),   # handles "2024-05-04"
                          parsed_time),
    parsed_time = if_else(is.na(parsed_time) & grepl("^\\d{4}$", timepoint),
                          ymd(paste0(timepoint, "-06-30")),   # assume midyear for year-only
                          parsed_time)
  ) %>%
  arrange(pool, parsed_time)


pool_time_summary <- subset(pool_time_summary, timepoint != "2023" & timepoint != "2024" & pool != "Joanna_coll")

# --- 2) Plot using parsed_time on x-axis ---
p <- ggplot(pool_time_summary,
            aes(x = parsed_time, y = mean_het, group = pool, color = pool)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_het - se_het, ymax = mean_het + se_het),
                width = 20) +   # width in days
  scale_x_date(date_labels = "%b %Y", date_breaks = "3 months") +
  labs(
    x = "Sampling date",
    y = "Mean heterozygosity (per-sample mean)",
    color = "Pool",
    title = "Average heterozygosity per pool over time"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )

ggsave("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/het_by_pool_timepoint_sorted.png",
       plot = p, width = 12, height = 8, dpi = 300)

















suppressPackageStartupMessages({
  library(dplyr)
  library(lubridate)
  library(ggplot2)
    library(dplyr)

})

# --- Per-sample stats (carry Group & QC fields so we can pick representatives)
het_per_sample <- het_merged %>%
  group_by(sample) %>%
  summarise(
    mean_het     = mean(het_prop, na.rm = TRUE),
    median_het   = median(het_prop, na.rm = TRUE),
    n_windows    = sum(!is.na(het_prop)),
    mean_depth   = mean(meanDepth, na.rm = TRUE),
    mean_missing = mean(missingRate, na.rm = TRUE),
    Group        = suppressWarnings(first(na.omit(Group))),  # one Group per sample
    .groups = "drop"
  )

# --- Merge metadata and build pool + parsed timepoint
md <- metadata_with_clone %>%
  filter(!(clone %in% c("Blank","BLANK"))) %>%
  mutate(pool_raw = ifelse(!is.na(accuratelocation) & nzchar(accuratelocation),
                           accuratelocation, location))

parse_to_date <- function(x) {
  d <- suppressWarnings(mdy(x))
  d <- ifelse(is.na(d), suppressWarnings(ymd(x)), d)
  d <- as.Date(d, origin = "1970-01-01")
  # year-only -> mid-year
  only_year <- grepl("^\\d{4}$", x)
  d[is.na(d) & only_year] <- ymd(paste0(x[is.na(d) & only_year], "-06-30"))
  as.Date(d)
}

merged <- het_per_sample %>%
  left_join(md, by = c("sample" = "Well")) %>%
  mutate(
    pool        = trimws(pool_raw),
    timepoint   = as.character(date),
    parsed_time = parse_to_date(timepoint)
  )
merged_slim <- merged %>%
  as_tibble() %>%
  transmute(
    pool         = as.character(pool),
    parsed_time  = as.Date(as.character(parsed_time)),
    Group        = as.character(Group),
    sample       = as.character(sample),
    n_windows    = as.numeric(n_windows),
    mean_depth   = as.numeric(mean_depth),
    mean_missing = as.numeric(mean_missing),
    mean_het     = as.numeric(mean_het),
    timepoint    = as.character(timepoint)
  ) %>%
  filter(!is.na(pool), !is.na(parsed_time), !is.na(Group))

# Pick ONE per (pool, parsed_time, Group)
# Tiebreakers: n_windows ↓, mean_depth ↓, mean_missing ↑, sample ↑
rep_one_per_group <- merged_slim %>%
  arrange(
    pool, parsed_time, Group,
    desc(n_windows), desc(mean_depth), mean_missing, sample
  ) %>%
  distinct(pool, parsed_time, Group, .keep_all = TRUE)

# (Optional) drop specific pools/timepoints
rep_one_per_group <- rep_one_per_group %>%
  filter(timepoint != "2023", timepoint != "2024", pool != "Joanna_coll")

# Aggregate PER POOL × TIME using the representatives only
pool_time_summary <- rep_one_per_group %>%
  group_by(pool, parsed_time) %>%
  summarise(
    n_groups   = n(),
    mean_het   = mean(mean_het, na.rm = TRUE),
    sd_het     = ifelse(n_groups > 1, sd(mean_het, na.rm = TRUE), 0),
    se_het     = ifelse(n_groups > 1, sd_het / sqrt(n_groups), NA_real_),
    median_het = median(mean_het, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(pool, parsed_time)


# --- Plot
p <- ggplot(pool_time_summary,
            aes(x = parsed_time, y = mean_het, group = pool, color = pool)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_het - se_het, ymax = mean_het + se_het), width = 20) +
  scale_x_date(date_labels = "%b %Y", date_breaks = "3 months") +
  labs(
    x = "Sampling date",
    y = "Mean heterozygosity (1 per Group per pool/timepoint)",
    color = "Pool",
    title = "Average heterozygosity per pool over time (group-representative samples)"
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank())

ggsave("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/het_by_pool_timepoint_grouprep.png",
       plot = p, width = 12, height = 8, dpi = 300)










suppressPackageStartupMessages({
  library(dplyr)
  library(lubridate)
  library(ggplot2)
    library(dplyr)

})

# --- Per-sample stats (carry Group & QC fields so we can pick representatives)

het_merged2 <- subset(het_merged, meanDepth >= 6)


het_per_sample <- het_merged2 %>%
  group_by(sample) %>%
  summarise(
    mean_het     = mean(het_prop, na.rm = TRUE),
    median_het   = median(het_prop, na.rm = TRUE),
    n_windows    = sum(!is.na(het_prop)),
    mean_depth   = mean(meanDepth, na.rm = TRUE),
    mean_missing = mean(missingRate, na.rm = TRUE),
    Group        = suppressWarnings(first(na.omit(Group))),  # one Group per sample
    .groups = "drop"
  )

# --- Merge metadata and build pool + parsed timepoint
md <- metadata_with_clone %>%
  filter(!(clone %in% c("Blank","BLANK"))) %>%
  mutate(pool_raw = ifelse(!is.na(accuratelocation) & nzchar(accuratelocation),
                           accuratelocation, location))

parse_to_date <- function(x) {
  d <- suppressWarnings(mdy(x))
  d <- ifelse(is.na(d), suppressWarnings(ymd(x)), d)
  d <- as.Date(d, origin = "1970-01-01")
  # year-only -> mid-year
  only_year <- grepl("^\\d{4}$", x)
  d[is.na(d) & only_year] <- ymd(paste0(x[is.na(d) & only_year], "-06-30"))
  as.Date(d)
}

merged <- het_per_sample %>%
  left_join(md, by = c("sample" = "Well")) %>%
  mutate(
    pool        = trimws(pool_raw),
    timepoint   = as.character(date),
    parsed_time = parse_to_date(timepoint)
  )
merged_slim <- merged %>%
  as_tibble() %>%
  transmute(
    pool         = as.character(pool),
    parsed_time  = as.Date(as.character(parsed_time)),
    Group        = as.character(Group),
    sample       = as.character(sample),
    n_windows    = as.numeric(n_windows),
    mean_depth   = as.numeric(mean_depth),
    mean_missing = as.numeric(mean_missing),
    mean_het     = as.numeric(mean_het),
    timepoint    = as.character(timepoint)
  ) %>%
  filter(!is.na(pool), !is.na(parsed_time), !is.na(Group))

# Pick ONE per (pool, parsed_time, Group)
# Tiebreakers: n_windows ↓, mean_depth ↓, mean_missing ↑, sample ↑
rep_one_per_group <- merged_slim %>%
  arrange(
    pool, parsed_time, Group,
    desc(n_windows), desc(mean_depth), mean_missing, sample
  ) %>%
  distinct(pool, parsed_time, Group, .keep_all = TRUE)

# (Optional) drop specific pools/timepoints
rep_one_per_group <- rep_one_per_group %>%
  filter(timepoint != "2023", timepoint != "2024", pool != "Joanna_coll")

# Aggregate PER POOL × TIME using the representatives only
pool_time_summary <- rep_one_per_group %>%
  group_by(pool, parsed_time) %>%
  summarise(
    n_groups   = n(),
    mean_het   = mean(mean_het, na.rm = TRUE),
    sd_het     = ifelse(n_groups > 1, sd(mean_het, na.rm = TRUE), 0),
    se_het     = ifelse(n_groups > 1, sd_het / sqrt(n_groups), NA_real_),
    median_het = median(mean_het, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(pool, parsed_time)


# --- Plot
p <- ggplot(pool_time_summary,
            aes(x = parsed_time, y = mean_het, group = pool, color = pool)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_het - se_het, ymax = mean_het + se_het), width = 20) +
  scale_x_date(date_labels = "%b %Y", date_breaks = "3 months") +
  labs(
    x = "Sampling date",
    y = "Mean heterozygosity (1 per Group per pool/timepoint)",
    color = "Pool",
    title = "Average heterozygosity per pool over time (group-representative samples)"
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank())

ggsave("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/het_by_pool_timepoint_grouprep_6andup.png",
       plot = p, width = 12, height = 8, dpi = 300)
