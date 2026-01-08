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
het_data <- read.table("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/Gilmer_het_100kb.txt",
                       header = FALSE, stringsAsFactors = FALSE)

colnames(het_data) <- c("contig", "win_start", "win_end", "sample", "het_prop")

# List all files











het_merged <- het_data %>%
  left_join(genomic_types, by = c("sample" = "CloneA")) %>%
  left_join(depths, by = c("sample" = "sampleId"))


het_mergedcont3 <- subset(het_merged, contig == "JAACYE010000003.1" & sample == "Rockpool4_D8")
het_mergedcont3_005 <- subset(het_mergedcont3, het_prop >= 0.005)

het_mergedcont3 <- subset(het_merged, contig == "JAACYE010000008.1" & sample == "Gilmer5_D6")
het_mergedcont3_005 <- subset(het_mergedcont3, het_prop >= 0.005)


het_mergedcont3 <- subset(het_merged, contig == "JAACYE010000011.1" & sample == "Gilmer5_D8")
het_mergedcont3_005 <- subset(het_mergedcont3, het_prop >= 0.005)


het_mergedcont3 <- subset(het_merged, contig == "JAACYE010000012.1" & sample == "Gilmer5_C9")
het_mergedcont3_005 <- subset(het_mergedcont3, het_prop >= 0.01)







"
contig win_start win_end
JAACYE010000003.1   2500000 2599999
JAACYE010000003.1   2600000 2699999

JAACYE010000007.1   1100000 1199999
JAACYE010000007.1   1200000 1299999

JAACYE010000008.1   1800000 1899999
JAACYE010000008.1   1900000 1999999

JAACYE010000011.1   1400000 1499999
JAACYE010000011.1   1500000 1599999
"



include_regions <- data.frame(
  contig = c(
    "JAACYE010000003.1", "JAACYE010000003.1",
 #   "JAACYE010000006.1", "JAACYE010000006.1",
    "JAACYE010000007.1", "JAACYE010000007.1",
    "JAACYE010000008.1", "JAACYE010000008.1",
    "JAACYE010000011.1", "JAACYE010000011.1",
    "JAACYE010000012.1"
  ),

  #need the bit on chr 6  Rockpool4_F7
  win_start = c(2500000, 2600000, 1100000, 1200000, 1800000, 1900000, 1400000, 1500000, 1100000),
  win_end   = c(2599999, 2699999, 1199999, 1299999, 1899999, 1999999, 1499999, 1599999, 1199999)
)

include_regions <- include_regions %>%
  mutate(
    contig = as.character(contig),
    win_start = as.integer(win_start),
    win_end   = as.integer(win_end)
  )

het_filtered_subset <- het_merged %>%
  mutate(
    contig = as.character(contig),
    win_start = as.integer(win_start),
    win_end   = as.integer(win_end)
  ) %>%
  # keep only exact matches of (contig, win_start, win_end)
  semi_join(include_regions, by = c("contig", "win_start", "win_end"))

# Check result
nrow(het_merged) - nrow(het_filtered_subset)

het_filtered_subset_3 <- subset(het_filtered_subset, meanDepth >= 4)


het_filtered_subset_3$norm_het <- het_filtered_subset_3$het_prop / sqrt(het_filtered_subset_3$meanDepth)


facet_means_fixed <- het_filtered_subset_3 %>%
  group_by(win_start, contig) %>%
  summarise(mean_norm_het = mean(het_prop, na.rm = TRUE))

subset(het_filtered_subset_3, contig == "JAACYE010000003.1")

het_filtered_subset_3 <- het_filtered_subset_3 %>%
  mutate(contigwindow = paste(contig, win_start, sep = "_"))

facet_means_fixed <- facet_means_fixed %>%
  mutate(contigwindow = paste(contig, win_start, sep = "_"))

facet_means_fixed$mean_het_fixed <- facet_means_fixed$mean_norm_het / mean(het_filtered_subset_3$meanDepth)

pdf("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/norm_het.pdf",
    width = 15, height = 10)

ggplot(het_filtered_subset_3, aes(x = norm_het, y = het_prop, color = meanDepth)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_vline(
    data = facet_means_fixed,
    aes(xintercept = mean_het_fixed),
    color = "red",
    linetype = "dashed",
    linewidth = 1,
    inherit.aes = FALSE
  ) +
  labs(
    x = "Normalized heterozygosity (norm_het)",
    y = "Raw heterozygosity (het_prop)",
    title = "Relationship between het_prop and norm_het across samples"
  ) +
  facet_wrap(~contigwindow, scales = "free_x", drop = TRUE) +
  theme_bw() +
  theme(text = element_text(size = 14))

dev.off()


library(dplyr)
library(tidyr)

df <- het_filtered_subset_3 %>%
  left_join(facet_means_fixed, by = "win_start")

# 2) Binary flag: 1 if above mean, 0 otherwise
df <- df %>%
  mutate(above_mean = if_else(norm_het > mean_het_fixed, 1L, 0L))

bin_mat <- df %>%
  select(win_start, sample, above_mean) %>%
  # If duplicates exist per (win_start, sample), collapse them to a single 0/1
  pivot_wider(
    names_from  = sample,
    values_from = above_mean,
    values_fill = list(above_mean = 0),
    values_fn   = list(above_mean = max)   # or ~ as.integer(any(. == 1))
  ) %>%
  arrange(win_start) %>%
  select(-win_start) %>%
  as.matrix()


signatures <- apply(bin_mat, 2, paste0, collapse = "")
sample_sign <- tibble(sample = colnames(bin_mat), signature = signatures)

# Count how many samples share each signature
sig_counts <- sample_sign %>%
  count(signature, name = "group_size")

# ---- 2) Label only those signatures that appear in >= 2 samples ----
# Excel-like A, B, ... Z, AA, AB, ...
num_to_letters <- function(n) {
  letter_index <- ((n - 1) %/% 26) + 1
  number_index <- ((n - 1) %% 26) + 1
  paste0(LETTERS[letter_index], "_", number_index)
}

multi_sig <- sig_counts %>%
  filter(group_size >= 2) %>%
  arrange(signature) %>%             # deterministic order
  mutate(label = sapply(seq_len(n()), num_to_letters))

# Map labels back to samples (singletons get NA label)
sample_labeled <- sample_sign %>%
  left_join(sig_counts, by = "signature") %>%
  left_join(select(multi_sig, signature, label), by = "signature") %>%
  select(sample, label, group_size, signature)

# ---- 3) Nice group summary (only labeled groups) ----
group_summary <- sample_labeled %>%
  filter(!is.na(label)) %>%
  group_by(label) %>%
  summarise(
    n = n(),
    members = paste(sample, collapse = ", "),
    .groups = "drop"
  ) %>%
  arrange(label)

# View results
print(sample_labeled)  # per-sample label (NA = no exact match with others)

sample_labeled <- as.data.frame(sample_labeled)

sample_labeled <- sample_labeled %>%
  mutate(label = ifelse(is.na(label), "A_XX", label))

group_summary <- sample_labeled %>%
  filter(!is.na(label)) %>%
  group_by(label) %>%
  summarise(
    n = n(),
    members = paste(sample, collapse = ", "),
    .groups = "drop"
  ) %>%
  arrange(label)

print(group_summary)   












pdf("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/norm_het2.pdf",
    width = 15, height = 10)

ggplot(het_filtered_subset_3, aes(x = contigwindow, y = het_prop, color = sample, group=sample)) +
  geom_point(size = 2, alpha = 0.7) + geom_line()+
  labs(
    x = "Contig and window",
    y = "Normalized heterozygosity (het_prop)",
    title = "Normalized het across windows (red = window mean)"
  ) + theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
 

dev.off()



het_filtered_subset_3label <- het_filtered_subset_3 %>%
  left_join(sample_labeled, by = c("sample")) 


pdf("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/norm_het3.pdf",
    width = 15, height = 10)
ggplot(het_filtered_subset_3label, aes(x = contigwindow, y = het_prop, color = label, group = sample)) +
  geom_line(alpha = 0.7) +
  geom_point(size = 2) +
  # Red mean trend line per contig
  geom_line(
    data = facet_means_fixed,
    aes(x = contigwindow, y = mean_het_fixed, group = contig),
    color = "red",
    linewidth = 1.2,
    inherit.aes = FALSE
  ) +
  # Red mean points
  geom_point(
    data = facet_means_fixed,
    aes(x = contigwindow, y = mean_het_fixed),
    color = "red",
    size = 3,
    inherit.aes = FALSE
  ) +
  labs(
    x = "Contig and window",
    y = "Normalized heterozygosity (het_prop)",
    title = "Normalized het across windows (red = window mean)"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    text = element_text(size = 14)
  )

dev.off()



pdf("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/norm_het_no_oo.pdf",
    width = 15, height = 10)

ggplot(subset(het_filtered_subset_3label, label != "A_XX"), aes(x = contigwindow, y = het_prop, color = label, group=sample)) +
  geom_point(size = 2, alpha = 0.7) + geom_line()+
  labs(
    x = "Contig and window",
    y = "Normalized heterozygosity (het_prop)",
    title = "Normalized het across windows (red = window mean)"
  )  + theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    text = element_text(size = 14)
  )

dev.off()




pdf("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/norm_het_only_oo.pdf",
    width = 15, height = 10)

ggplot(subset(het_filtered_subset_3label, label == "A_XX"), aes(x = contigwindow, y = het_prop, color = sample, group=sample)) +
  geom_point(size = 2, alpha = 0.7) + geom_line()+
  labs(
    x = "Contig and window",
    y = "Normalized heterozygosity (het_prop)",
    title = "Normalized het across windows (red = window mean)"
  )  + theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    text = element_text(size = 14)
  )

dev.off()


metadata_with_cloneA_label <- metadata_with_clone %>%
  left_join(sample_labeled, by = c("Well" = "sample")) 





metadata_with_cloneA_label <- subset(metadata_with_cloneA_label,
                                     !(date %in% c("2023", "2024")))
metadata_with_cloneA_label_no_NA <- metadata_with_cloneA_label %>%
  filter(!is.na(group_size))


metadata_with_cloneA_label_no_NA_filtered <- metadata_with_cloneA_label_no_NA %>%
  select(Well, date, label, accuratelocation, signature)



# Count (sum) by superclone and date
subgroups_count <- metadata_with_cloneA_label_no_NA_filtered %>%
  group_by(label, accuratelocation, date, signature) %>%
  summarise(count = n(), .groups = "drop")


subgroups_count <- as.data.frame(subgroups_count)

subgroups_count$date <- as.Date(subgroups_count$date, format = "%m/%d/%Y")
subgroups_count$date[is.na(subgroups_count$date)] <- as.Date(paste0(subgroups_count$date[is.na(subgroups_count$date)], "-01-01"))

# Reorder factor levels by date
subgroups_count$date <- factor(subgroups_count$date, levels = unique(subgroups_count$date[order(subgroups_count$date)]))


library(ggalluvial)



png("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/Gilmer_counts_signature.png", res = 600,
    width = 8000, height = 4000)

# Plot
ggplot(subgroups_count,
       aes(x = date, 
           stratum = signature, 
           alluvium = signature, 
           y = count,
           fill = signature,
           label = signature)) +
  geom_flow(alpha = 0.5, width = 0.3) +
  geom_stratum(width = 0.3, color = "grey30") +
  geom_text(stat = "stratum", size = 3) +
  scale_fill_viridis_d(option = "C") +
  theme_bw() + facet_grid(~accuratelocation, scales = "free_x")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Gilmer subtypes Transitions Over Time",
       y = "Count",
       x = "Timepoint")

dev.off()





subgroups_plot <- subgroups_count %>%
  mutate(date = as.Date(date)) %>%                 # ensure Date type
  group_by(accuratelocation, date, label) %>%
  summarise(count = sum(count), .groups = "drop") %>%
  arrange(accuratelocation, date)

# 2) Make a per-facet ordered x-axis (factor), so ggalluvial sees discrete axes
subgroups_plot <- subgroups_plot %>%
  group_by(accuratelocation) %>%
  mutate(date_fac = factor(format(date), levels = unique(format(date)))) %>%
  ungroup()

png("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/Gilmer_counts.png",
    res = 600, width = 8000, height = 4000)

ggplot(
  subgroups_plot,
  aes(x = date_fac,
      stratum = label,        # strata shown at each timepoint
      alluvium = label,       # same id across time so flows connect
      y = count,
      fill = label,
      label = label)
) +
  geom_flow(alpha = 0.5, width = 0.3) +
  geom_stratum(width = 0.3, color = "grey30") +
  geom_text(stat = "stratum", size = 3) +
  scale_fill_viridis_d(option = "C") +
  facet_grid(~ accuratelocation, scales = "free_x", drop = TRUE) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Gilmer subtypes transitions over time",
    y = "Count",
    x = "Timepoint"
  )

dev.off()






















# Plot
pdf("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/het_window_Gilmer_new_100000.pdf",
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







  p2 <- ggplot(het_merged, aes(x = win_start, y = het_prop, group= sample, col=meanDepth)) +
    geom_line() +
    facet_grid(~contig, scales = "free_x") +
    labs(
      title = "Windowed heterozygosity (100kb)",
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
pdf("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/het_window_Gilmer_100000_overlap.pdf",
    width = 100, height = 10)
p2
dev.off()




  p3 <- ggplot(subset(het_merged, meanDepth > 5), aes(x = win_start, y = het_prop, group= sample, col=meanDepth)) +
    geom_line() +
    facet_grid(~contig, scales = "free_x") +
    labs(
      title = "Windowed heterozygosity (100kb)",
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
pdf("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/het_window_Gilmer_100000_overlap_depth5.pdf",
    width = 100, height = 10)
p3
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
pdf("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/het_window_Gilmer_100000_overlap_depth3.pdf",
    width = 100, height = 10)
p3
dev.off()







het_data <- read.table("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/Gilmer_het_100kb.txt",
                       header = FALSE, stringsAsFactors = FALSE)

colnames(het_data) <- c("contig", "win_start", "win_end", "sample", "het_prop")

# List all files
het_data <- subset(het_data, contig == "JAACYE010000003.1")

het_merged <- het_data %>%
  left_join(genomic_types, by = c("sample" = "CloneA")) %>%
  left_join(depths, by = c("sample" = "sampleId")) 


het_mergedgilmer <- subset(het_merged, meanDepth >= 3)

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

png("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/pca_Gilmer_A_100kb_het.png",
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



png("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/pca_Gilmer_A_100kb_het_2_3.png",
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




png("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/pca_Gilmer_A_100kb_het_3_4.png",
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
