#ijob -A berglandlab -c10 -p standard --mem=40G
#module load gcc/11.4.0  openmpi/4.1.4 icu R/4.3.1
#R

library(data.table)
fst_raw <- read.table("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/results_pixy10000/pixy_10000__fst.txt",header = TRUE,sep = "\t")
dxy_raw <- read.table("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/results_pixy10000/pixy_10000__dxy.txt",header = TRUE,sep = "\t")
pi_raw <- read.table("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/results_pixy10000/pixy_10000__pi.txt",header = TRUE,sep = "\t")

library(ggplot2)


chrom_keep <- paste0("JAACYE0100000", sprintf("%02d", 1:12), ".1")

# Subset the data
fst_subset <- fst_raw[fst_raw$chromosome %in% chrom_keep, ]

fst_subsetsubset <- subset(fst_subset, chromosome == "JAACYE010000007.1")

fst_1060001to1120001 <- subset(fst_subsetsubset, window_pos_1 > 1060000 & window_pos_1 < 1120002)
unique(fst_subsetsubset$window_pos_1)
# Check
unique(fst_subset$chromosome)

# Clean up
fst_clean <- na.omit(fst_subset)


fst_clean$window_pos_1 <- as.numeric(fst_clean$window_pos_1)
fst_clean$avg_wc_fst   <- as.numeric(fst_clean$avg_wc_fst)

# Add combined comparison label
fst_clean$comparison <- paste(fst_clean$pop1, fst_clean$pop2, sep = "_")

# Save plot as PNG
png("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/fst_manhattan10000.png",
    width = 10000, height = 5000, res = 300)

ggplot(fst_clean, aes(x = window_pos_1, y = avg_wc_fst)) +
  geom_point(size = 0.5, alpha = 0.6, color = "darkblue") +
  facet_grid(comparison ~ chromosome, scales = "free_x", space = "free_x") +
  labs(
    title = "Manhattan-style FST Plot",
    x = "Genomic position (window start)",
    y = "FST"
  ) + ylim(0,1)+
  theme_bw(base_size = 12) +
  theme(
    strip.text.x = element_text(size = 10, angle = 90),
    strip.text.y = element_text(size = 8),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.minor = element_blank()
  )

dev.off()








chrom_keep <- paste0("JAACYE0100000", sprintf("%02d", 1:12), ".1")

dxy_subset <- dxy_raw[dxy_raw$chromosome %in% chrom_keep, ]


# Check
unique(dxy_subset$chromosome)

# Clean up
dxy_clean <- na.omit(dxy_subset)


dxy_clean$window_pos_1 <- as.numeric(dxy_clean$window_pos_1)
dxy_clean$avg_dxy   <- as.numeric(dxy_clean$avg_dxy)

# Add combined comparison label
dxy_clean$comparison <- paste(dxy_clean$pop1, dxy_clean$pop2, sep = "_")

# Save plot as PNG
png("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/dxy_manhattan10000.png",
    width = 10000, height = 5000, res = 300)

ggplot(dxy_clean, aes(x = window_pos_1, y = avg_dxy)) +
  geom_point(size = 0.5, alpha = 0.6, color = "darkblue") +
  facet_grid(comparison ~ chromosome, scales = "free_x", space = "free_x") +
  labs(
    title = "Manhattan-style dxy Plot",
    x = "Genomic position (window start)",
    y = "dxy"
  ) + ylim(0,0.1)+
  theme_bw(base_size = 12) +
  theme(
    strip.text.x = element_text(size = 10, angle = 90),
    strip.text.y = element_text(size = 8),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.minor = element_blank()
  )

dev.off()






















chrom_keep <- paste0("JAACYE0100000", sprintf("%02d", 1:12), ".1")

pi_subset <- pi_raw[pi_raw$chromosome %in% chrom_keep, ]


# Check
unique(pi_subset$chromosome)

# Clean up
pi_clean <- na.omit(pi_subset)


pi_clean$window_pos_1 <- as.numeric(pi_clean$window_pos_1)
pi_clean$avg_pi   <- as.numeric(pi_clean$avg_pi)

# Add combined comparison label

# Save plot as PNG
png("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/pi_manhattan10000.png",
    width = 10000, height = 5000, res = 300)

ggplot(pi_clean, aes(x = window_pos_1, y = avg_pi)) +
  geom_point(size = 0.5, alpha = 0.6, color = "darkblue") +
  facet_grid(pop ~ chromosome, scales = "free_x", space = "free_x") +
  labs(
    title = "Manhattan-style pi Plot",
    x = "Genomic position (window start)",
    y = "pi"
  ) + ylim(0,0.1)+
  theme_bw(base_size = 12) +
  theme(
    strip.text.x = element_text(size = 10, angle = 90),
    strip.text.y = element_text(size = 8),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.minor = element_blank()
  )

dev.off()







png("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/dxy_box_10000.png",
    width = 5000, height = 2000, res = 300)

ggplot(dxy_clean, aes(x = comparison, y = avg_dxy)) +
  geom_boxplot() +
  labs(
    title = "Manhattan-style dxy Plot",
    x = "Genomic position (window start)",
    y = "dxy"
  ) + ylim(0,0.1)+
  theme_bw(base_size = 12) +
  theme(
    strip.text.x = element_text(size = 10, angle = 90),
    strip.text.y = element_text(size = 8),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

dev.off()




png("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/pi_box_10000.png",
    width = 5000, height = 2000, res = 300)

ggplot(pi_clean, aes(x = pop, y = avg_pi)) +
  geom_boxplot() +
  labs(
    title = "pi box plot Plot",
    x = "pop",
    y = "pi"
  ) + ylim(0,0.1)+
  theme_bw(base_size = 12) +
  theme(
    strip.text.x = element_text(size = 10, angle = 90),
    strip.text.y = element_text(size = 8),
    panel.grid.minor = element_blank()
  )

dev.off()






library(dplyr)
library(purrr)
library(tibble)
library(readr)
library(tidyr)



library(dplyr)
library(purrr)
library(tidyr)




# Example: define your make_windows function
make_windows <- function(chr, len, win_size = 50000) {
  tibble(
    chromosome = chr,
    window_start = seq(1, len, by = win_size),
    window_end   = pmin(seq(1, len, by = win_size) + win_size - 1, len)
  )
}

# 1. Generate windows per chromosome
windows_all <- chr_lengths %>%
  mutate(windows = map2(chromosome, chr_len, make_windows)) %>%
  select(-chr_len) %>%
  unnest(windows, names_sep = "_") %>%   # prevent column name clashes
  rename(
    chromosome = windows_chromosome,
    window_start = windows_window_start,
    window_end   = windows_window_end
  )

# 2. Expand windows for each population
windows_expanded <- expand_grid(
  pop = unique(pi_clean$pop),
  windows_all
)

# 3. Summarise pi into windows
df_binned <- windows_expanded %>%
  rowwise() %>%   # safe to evaluate per row
  mutate(
    avg_pi = {
      sub <- pi_clean %>%
        filter(
          pop == pop,
          chromosome == chromosome,
          window_pos_1 >= window_start,
          window_pos_2 <= window_end
        )
      if (nrow(sub) == 0) NA_real_
      else weighted.mean(sub$avg_pi, sub$no_sites)
    },
    no_sites = {
      sub <- pi_clean %>%
        filter(
          pop == pop,
          chromosome == chromosome,
          window_pos_1 >= window_start,
          window_pos_2 <= window_end
        )
      sum(sub$no_sites, na.rm = TRUE)
    }
  ) %>%
  ungroup()


  df_binned2 <- as.data.frame(df_binned)

  head(df_binned2)  
  tail(df_binned2)





  
png("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/pi_box_10000_windowed_50000.png",
    width = 5000, height = 2000, res = 300)

ggplot(df_binned2, aes(x = pop, y = avg_pi)) +
  geom_boxplot() +
  labs(
    title = "pi box plot",
    x = "pop",
    y = "pi"
  ) + ylim(0,0.1)+
  theme_bw(base_size = 12) +
  theme(
    strip.text.x = element_text(size = 10, angle = 90),
    strip.text.y = element_text(size = 8),
    panel.grid.minor = element_blank()
  )

dev.off()








