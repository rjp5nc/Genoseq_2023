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



het_win <- het_merged %>%
  mutate(
    contig    = as.character(contig),
    win_start = as.integer(win_start),
    win_end   = as.integer(win_end)
  ) %>%
  filter(meanDepth >= 4)

# Optional normalization by depth (comment out if you do not want it)
het_win <- het_win %>%
  mutate(norm_het = het_prop / sqrt(meanDepth))

# Canonical window id and ordering
het_win <- het_win %>%
  mutate(window_id = paste(contig, sprintf("%010d", win_start), sep = "_"))

win_levels <- het_win %>%
  distinct(contig, win_start, win_end, window_id) %>%
  arrange(contig, win_start) %>%
  pull(window_id)

het_win <- het_win %>%
  mutate(window_id = factor(window_id, levels = win_levels))
window_means2 <- het_win %>%
  group_by(contig, win_start) %>%
  summarise(mean_het = mean(het_prop, na.rm = TRUE), .groups = "drop")

# Ensure plotting order within each group
het_win_plot <- het_win %>%
  arrange(contig, sample, win_start)

pdf("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/het_per_window.pdf",
    width = 15, height = 80)

ggplot(
  het_win_plot,
  aes(x = win_start, y = het_prop,
      color = sample, group = interaction(sample, contig))
) +
  geom_line(alpha = 0.6) +
  geom_point(size = 0.8) +
  geom_line(
    data = window_means2,
    aes(x = win_start, y = mean_het, group = 1),
    inherit.aes = FALSE, color = "red", linewidth = 0.9
  ) +
  geom_point(
    data = window_means2,
    aes(x = win_start, y = mean_het),
    inherit.aes = FALSE, color = "red", size = 1.8
  ) +
  facet_wrap(. ~ contig, scales = "free_x", ncol=1) +
  labs(
    x = "Window start (bp)",
    y = "het_prop",
    title = "Per-window heterozygosity by contig (red = per-window mean)"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    text = element_text(size = 14),
    panel.grid.major.x = element_blank()
  )
dev.off()

# --------------- Above/below mean flags at window level -----------------------
df_w <- het_win %>%
  left_join(window_means, by = c("contig", "win_start", "win_end", "window_id")) %>%
  mutate(above_mean = if_else(het_prop > mean_het, 1L, 0L))

# Binary matrix: rows = windows, cols = samples
bin_mat_w <- df_w %>%
  select(window_id, sample, above_mean) %>%
  distinct() %>%
  pivot_wider(names_from = sample, values_from = above_mean,
              values_fill = list(above_mean = 0)) %>%
  arrange(window_id) %>%
  tibble::column_to_rownames("window_id") %>%
  as.matrix()

# Jaccard similarity across samples
jaccard_sim <- function(M) {
  p <- ncol(M)
  S <- matrix(NA_real_, p, p, dimnames = list(colnames(M), colnames(M)))
  for (i in seq_len(p)) for (j in seq_len(p)) {
    a <- M[, i] > 0; b <- M[, j] > 0
    inter <- sum(a & b); uni <- sum(a | b)
    S[i, j] <- if (uni == 0) 1 else inter / uni
  }
  S
}
sim_matrix_w <- jaccard_sim(bin_mat_w)

# ------------------- Reuse your signature grouping (optional) -----------------
# Exact-match signature per sample using the per-window flags
signatures <- apply(bin_mat_w, 2, paste0, collapse = "")
sample_sign <- tibble(sample = colnames(bin_mat_w), signature = signatures)
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

# If you want the signature-colored plot (per window) --------------------------
het_win_sig <- het_win %>%
  left_join(sample_labeled, by = c("sample" = "sample"))

het_win_sig <- het_win_sig %>% filter(!is.na(signature))

pdf("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/het_per_window_signature.pdf",
    width = 15, height = 10)

ggplot(het_win_sig, aes(x = window_id, y = het_prop, color = signature, group = sample)) +
  geom_line(alpha = 0.6) +
  geom_point(size = 1.2) +
  geom_line(data = window_means,
            aes(x = window_id, y = mean_het, group = 1),
            inherit.aes = FALSE, color = "red", linewidth = 1.0) +
  geom_point(data = window_means,
             aes(x = window_id, y = mean_het),
             inherit.aes = FALSE, color = "red", size = 2) +
  labs(x = "Genome windows (by contig, by start)",
       y = "het_prop",
       title = "Per-window heterozygosity across the genome by signature") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        text = element_text(size = 14))

dev.off()
