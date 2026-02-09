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


het_mergedcont3 <- subset(het_merged, contig == "JAACYE010000003.1")
het_mergedcont3_005 <- subset(het_mergedcont3, het_prop >= 0.01)
het_mergedcont3_005 

het_mergedcont3 <- subset(het_merged, contig == "JAACYE010000008.1" & sample == "Gilmer5_D6")
het_mergedcont3_005 <- subset(het_mergedcont3, het_prop >= 0.005)


het_mergedcont3 <- subset(het_merged, contig == "JAACYE010000002.1")
het_mergedcont3_005 <- subset(het_mergedcont3, het_prop >= 0.01)
het_mergedcont3_005

het_mergedcont3 <- subset(het_merged, contig == "JAACYE010000012.1" & sample == "Gilmer5_C9")
het_mergedcont3_005 <- subset(het_mergedcont3, het_prop >= 0.01)




0.0015


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


window_blocks <- het_filtered_subset_3 %>%
  distinct(contig, win_start, win_end) %>%
  arrange(contig, win_start) %>%
  group_by(contig) %>%
  mutate(block = cumsum(is.na(lag(win_end)) | (win_start != lag(win_end) + 1))) %>%
  ungroup()

# Attach block id back to each row
het_with_blocks <- het_filtered_subset_3 %>%
  inner_join(window_blocks, by = c("contig", "win_start", "win_end")) %>%
  mutate(contig_block = paste0(contig, "_b", sprintf("%02d", block)))

# --- 2) Collapse ONLY within each contiguous block (per sample) ---
het_collapsed <- het_with_blocks %>%
  group_by(contig, block, contig_block, sample) %>%
  summarise(
    mean_het_prop = mean(het_prop, na.rm = TRUE),
    meanDepth_c   = mean(meanDepth, na.rm = TRUE),
    .groups = "drop"
  ) 

# Optional: order contig_block nicely (by contig, then block)
cb_levels <- het_collapsed %>%
  distinct(contig, block, contig_block) %>%
  arrange(contig, block) %>%
  pull(contig_block)

het_collapsed <- het_collapsed %>%
  mutate(contig_block = factor(contig_block, levels = cb_levels))

# --- 3) Per-block means (for thresholds / red line) ---

block_means <- het_collapsed %>%
  group_by(contig, block, contig_block) %>%
  summarise(
    mean_norm_het_c = 0.0015,
    .groups = "drop")


# --- 4) Plot (one value per contiguous block per sample) + red mean line/points ---
pdf("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/NOTNORM_het_contig_blocks.pdf",
    width = 15, height = 10)

ggplot(subset(het_collapsed), aes(x = contig_block, y = mean_het_prop, color = sample, group = sample)) +
  geom_line(alpha = 0.7) +
  geom_point(size = 2) +
  # red mean line per contig (segments between blocks of same contig)
  geom_line(
    data = block_means,
    aes(x = contig_block, y = mean_norm_het_c, group = contig),
    color = "red", linewidth = 1.2, inherit.aes = FALSE
  ) +
    geom_point(
    data = block_means,
    aes(x = contig_block, y = mean_norm_het_c),
    color = "red", size = 3, inherit.aes = FALSE
  ) +
  labs(
    x = "Contig (contiguous block)",
    y = "Mean het_prop (within contiguous window block)",
    title = "Collapsed het across contiguous blocks (red = per-block mean)"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        text = element_text(size = 14))

dev.off()




# --- 5) Above/below-mean flag (per contiguous block) ---
df_c <- het_collapsed %>%
  left_join(block_means, by = c("contig", "block", "contig_block")) %>%
  mutate(above_mean = if_else(mean_het_prop > mean_norm_het_c, 1L, 0L))

# --- 6) Similarity matrix over contiguous blocks (rows) × samples (cols) ---
bin_mat_c <- df_c %>%
  select(contig_block, sample, above_mean) %>%
  distinct() %>%
  pivot_wider(
    names_from  = sample,
    values_from = above_mean,
    values_fill = list(above_mean = 0)
  ) %>%
  arrange(contig_block) %>%
  tibble::column_to_rownames("contig_block") %>%
  as.matrix()

# Jaccard similarity across samples
jaccard_sim <- function(M) {
  p <- ncol(M)
  S <- matrix(NA_real_, p, p, dimnames = list(colnames(M), colnames(M)))
  for (i in seq_len(p)) for (j in seq_len(p)) {
    a <- M[, i]; b <- M[, j]
    inter <- sum(a & b); uni <- sum(a | b)
    S[i, j] <- if (uni == 0) 1 else inter / uni
  }
  S
}
sim_matrix_c <- jaccard_sim(bin_mat_c)

# (Optional) exact-match grouping on collapsed blocks (reuse your labeling)
signatures <- apply(bin_mat_c, 2, paste0, collapse = "")
sample_sign <- tibble(sample = colnames(bin_mat_c), signature = signatures)
sig_counts  <- sample_sign %>% count(signature, name = "group_size")

num_to_letters <- function(n) {
  letter_index <- ((n - 1) %/% 26) + 1
  number_index <- ((n - 1) %% 26) + 1
  paste0(LETTERS[letter_index], "_", number_index)
}
multi_sig <- sig_counts %>%
  filter(group_size >= 2) %>%
  arrange(signature) %>%
  mutate(label = sapply(seq_len(n()), num_to_letters))

sample_labeled <- sample_sign %>%
  left_join(sig_counts, by = "signature") %>%
  left_join(select(multi_sig, signature, label), by = "signature") %>%
  mutate(label = if_else(is.na(label), "A_XX", label)) %>%
  select(sample, label, group_size, signature)


sample_labeled <- as.data.frame(sample_labeled)




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








hetcollapsed_graph <- as.data.frame(het_collapsed) %>%
  left_join(metadata_with_cloneA_label_no_NA_filtered, by = c("sample" = "Well")) 


hetcollapsed_graph <- as.data.frame(hetcollapsed_graph)
hetcollapsed_graph <- hetcollapsed_graph %>% filter(!is.na(signature))

# --- 4) Plot (one value per contiguous block per sample) + red mean line/points ---
pdf("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/NOTnorm_het_contig_blocks_signature.pdf",
    width = 15, height = 10)

ggplot(hetcollapsed_graph, aes(x = contig_block, y = mean_het_prop, color = signature, group = sample)) +
  geom_line(alpha = 0.7) +
  geom_point(size = 2) +
  # red mean line per contig (segments between blocks of same contig)
  geom_line(
    data = block_means,
    aes(x = contig_block, y = mean_norm_het_c, group = contig),
    color = "red", linewidth = 1.2, inherit.aes = FALSE
  ) +
  geom_point(
    data = block_means,
    aes(x = contig_block, y = mean_norm_het_c),
    color = "red", size = 3, inherit.aes = FALSE
  ) +
  labs(
    x = "Contig (contiguous block)",
    y = "Mean het_prop (within contiguous window block)",
    title = "Collapsed het across contiguous blocks (red = per-block mean)"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        text = element_text(size = 14))

dev.off()


write.csv(hetcollapsed_graph, "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/Gilmer_A_het_blocks_signature.csv")
write.csv(metadata_with_cloneA_label_no_NA_filtered, "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/Gilmer_A_het_blocks_all.csv")





# Filter for signature 01111101
sig_target <- "01111101"
selected_individuals <- hetcollapsed_graph %>%
  filter(signature == sig_target) %>%
  distinct(sample)

# Save list of individuals to a text file
out_list <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/Gilmer_signature_01111101_samples.txt"
write.table(selected_individuals$sample, file = out_list,
            quote = FALSE, row.names = FALSE, col.names = FALSE)

cat("Saved sample list to:", out_list, "\n")





subgroups_clean <- subgroups_count %>%
  mutate(
    signature = sprintf("%08s", signature),
    signature = gsub(" ", "0", signature),
    date      = as.Date(date)  # assumes YYYY-MM-DD or something parseable
  )

# 2) Split signature into bit1..bit8
sub_bits <- subgroups_clean %>%
  mutate(sig_split = strsplit(signature, "")) %>%
  unnest_wider(sig_split, names_sep = "") %>%
  rename(
    bit1 = sig_split1,
    bit2 = sig_split2,
    bit3 = sig_split3,
    bit4 = sig_split4,
    bit5 = sig_split5,
    bit6 = sig_split6,
    bit7 = sig_split7,
    bit8 = sig_split8
  )

bit_cols <- paste0("bit", 1:8)

# 3) Long format: one row per (date, bit_position, bit_value)
bit_long <- sub_bits %>%
  pivot_longer(
    cols      = all_of(bit_cols),
    names_to  = "bit_position",
    values_to = "bit_value"
  ) %>%
  mutate(
    bit_position = as.numeric(gsub("bit", "", bit_position)),
    bit_value    = as.numeric(bit_value)
  )

# 4) Aggregate counts and compute percentages
bit_pct <- bit_long %>%
  group_by(date, bit_position, bit_value) %>%
  summarise(total_count = sum(count), .groups = "drop") %>%
  group_by(date, bit_position) %>%
  mutate(percentage = 100 * total_count / sum(total_count)) %>%
  ungroup()

# 5) Ensure both bit_value 0 and 1 exist for each (date, bit_position)
bit_pct_full <- bit_pct %>%
  mutate(bit_value = factor(bit_value, levels = c(0, 1))) %>%
  complete(
    date,
    bit_position,
    bit_value,
    fill = list(total_count = 0, percentage = 0)
  )


  bit_totals <- bit_pct %>%
  group_by(bit_position, bit_value) %>%       # ignore date
  summarise(
    total_count = sum(total_count),
    .groups = "drop"
  ) %>%
  group_by(bit_position) %>%
  mutate(
    pct = 100 * total_count / sum(total_count)
  ) %>%
  ungroup()


bit_totals$bit_value <- as.character(bit_totals$bit_value)

  bitsovertime_time_total <- ggplot(
  bit_totals,
  aes(
    x     = bit_position,
    y     = total_count,
    color = bit_value,
    group = interaction(bit_value, bit_position)
  )
) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2)
  labs(
    title = "Percentage of 0's and 1's per bit over time",
    x     = "Sampling date",
    y     = "Percentage (%)"
  ) +
  theme_bw() +
  theme(
    strip.text  = element_text(size = 12),
    plot.title  = element_text(face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )

ggsave(
  "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/outbits_time_totalcount.png",
  bitsovertime_time_total,
  width  = 6,
  height = 4,
  dpi    = 300
)







  bitsovertime_time <- ggplot(
  bit_pct_full,
  aes(
    x     = date,
    y     = percentage,
    color = bit_value,
    group = interaction(bit_value, bit_position)
  )
) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2) +
  scale_color_manual(
    values = c("0" = "dodgerblue3", "1" = "firebrick2"),
    name   = "Bit value"
  ) +
  scale_x_date(date_labels = "%Y-%m-%d") +
  facet_wrap(~ bit_position, ncol = 4) +
  labs(
    title = "Percentage of 0's and 1's per bit over time",
    x     = "Sampling date",
    y     = "Percentage (%)"
  ) +
  theme_bw() +
  theme(
    strip.text  = element_text(size = 12),
    plot.title  = element_text(face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )

ggsave(
  "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/outbits_time.png",
  bitsovertime_time,
  width  = 12,
  height = 8,
  dpi    = 300
)







bit_pct_het <- subset(bit_pct_full, bit_value == 1)

bit_pct_het$bit_position <- as.character(bit_pct_het$bit_position)

  bitsovertime_het <- ggplot(
  bit_pct_het,
  aes(
    x     = date,
    y     = percentage,
    color = bit_position,
    group = bit_position
  )
) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2) +
  scale_x_date(date_labels = "%Y-%m-%d") +
  labs(
    title = "Percentage of 1's per bit over time (heterozygosity over time)",
    x     = "Sampling date",
    y     = "Percentage (%)"
  ) +
  theme_bw() 

ggsave(
  "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/outbits_time_together.png",
  bitsovertime_het,
  width  = 12,
  height = 8,
  dpi    = 300
)






bitsovertime_bydate <- ggplot(
  bit_pct_full,
  aes(
    x     = bit_position,
    y     = percentage,
    color = bit_value,
    group = bit_value
  )
) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = 1:8) +
  scale_color_manual(
    values = c("0" = "dodgerblue3", "1" = "firebrick2"),
    name   = "Bit value"
  ) +
  facet_wrap(~ date, ncol = 1) +
  labs(
    title = "Percentage of 0's and 1's across bit positions",
    x     = "Bit position (1-8)",   # plain hyphen to avoid encoding warning
    y     = "Percentage (%)"
  ) +
  theme_bw() +
  theme(
    strip.text  = element_text(size = 12),
    plot.title  = element_text(face = "bold", size = 14),
    legend.position = "top"
  )

ggsave(
  "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/outbits_bydate.png",
  bitsovertime_bydate,
  width  = 10,
  height = 12,
  dpi    = 300
)

