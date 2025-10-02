#ijob -A berglandlab -c10 -p standard --mem=40G

#module load gcc/11.4.0  openmpi/4.1.4 icu R/4.3.1
#R

library(tidyverse)
library(ggplot2)
library(data.table)
library(patchwork)

dp_raw <- read.table("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/avgdepth_all_100000.txt",
                     header = FALSE, sep = " ")

colnames(dp_raw) <- c("contig", "window_start", "window_end", "sample", "depth")

dp_clean <- na.omit(dp_raw)


dp_clean$window_start <- as.numeric(dp_clean$window_start)
dp_clean$depth   <- as.numeric(dp_clean$depth)


# Save plot as PNG
png("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/dp_manhattan_100k.png",
    width = 10000, height = 1000, res = 300)

ggplot(dp_clean, aes(x = window_start, y = avg_depth)) +
  geom_point(size = 0.5, alpha = 0.6, color = "darkblue") +
  facet_grid(~contig, scales = "free_x", space = "free_x") +
  labs(
    title = "Manhattan-style dp Plot",
    x = "Genomic position (window start)",
    y = "dp"
  ) + ylim(0,20)+
  theme_bw(base_size = 12) +
  theme(
    strip.text.x = element_text(size = 10, angle = 90),
    strip.text.y = element_text(size = 8),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.minor = element_blank()
  )

dev.off()



# Save plot as PNG
png("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/dp_7_manhattan_100k.png",
    width = 10000, height = 1000, res = 300)

ggplot(subset(dp_clean, contig == "JAACYE010000007.1"), aes(x = window_start, y = avg_depth)) +
  geom_point(size = 0.5, alpha = 0.6, color = "darkblue") +
  facet_grid(~contig) +
  labs(
    title = "Manhattan-style dp Plot",
    x = "Genomic position (window start)",
    y = "dp"
  ) + ylim(0,20)+
    scale_x_continuous(breaks = seq(0, max(dp_clean$window_start), by = 1000000)) +
  theme_bw(base_size = 12) 

dev.off()

















#dp_raw_all <- read.table("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/avg.depth.all12.txt",
#                     header = FALSE, sep = "\t")

dp_clean_all <- na.omit(dp_raw_all)

colnames(dp_raw_all) <- c("contig", "window_start", "window_end", "sample", "depth")


dp_clean_all$window_start <- as.numeric(dp_clean_all$window_start)
dp_clean_all$depth   <- as.numeric(dp_clean_all$depth)

pdf("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/dp_manhattan_by_sample_100k.pdf",
    width = 25, height = 2)

for (s in unique(dp_clean_all$sample)) {
  p <- ggplot(subset(dp_clean_all, sample == s),
              aes(x = window_start, y = depth)) +
    geom_point(size = 0.4, alpha = 0.6, color = "darkblue") +
    facet_grid(~contig, scales = "free_x") +
    labs(
      title = paste("Depth profile:", s),
      x = "Genomic position (window start)",
      y = "Depth"
    ) +
    ylim(0, 40) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.minor = element_blank()
    )
  print(p)
}

dev.off()


dp_clean_all <- na.omit(dp_raw)

colnames(dp_clean_all) <- c("contig", "window_start", "window_end", "sample", "depth")


dp_clean_all$window_start <- as.numeric(dp_clean_all$window_start)
dp_clean_all$depth   <- as.numeric(dp_clean_all$depth)

pdf("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/dp_manhattan_by_sample_100k.pdf",
    width = 25, height = 2)

for (s in unique(dp_clean_all$sample)) {
  p <- ggplot(subset(dp_clean_all, sample == s),
              aes(x = window_start, y = depth)) +
    geom_point(size = 0.4, alpha = 0.6, color = "darkblue") +
    facet_grid(~contig, scales = "free_x") +
    labs(
      title = paste("Depth profile:", s),
      x = "Genomic position (window start)",
      y = "Depth"
    ) +
    ylim(0, 40) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.minor = element_blank()
    )
  print(p)
}

dev.off()






dp_raw_all <- read.table("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/avgdepth_all_100000.txt",
                     header = FALSE, sep = " ")

dp_clean_all <- na.omit(dp_raw_all)

colnames(dp_clean_all) <- c("contig", "window_start", "window_end", "sample", "depth")


dp_clean_all$window_start <- as.numeric(dp_clean_all$window_start)
dp_clean_all$depth   <- as.numeric(dp_clean_all$depth)

pdf("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/dp_manhattan_by_sample_100k.pdf",
    width = 25, height = 2)

for (s in unique(dp_clean_all$sample)) {
  p <- ggplot(subset(dp_clean_all, sample == s),
              aes(x = window_start, y = depth)) +
    geom_point(size = 0.4, alpha = 0.6, color = "darkblue") +
    facet_grid(~contig, scales = "free_x") +
    labs(
      title = paste("Depth profile:", s),
      x = "Genomic position (window start)",
      y = "Depth"
    ) +
    ylim(0, 40) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.minor = element_blank()
    )
  print(p)
}

dev.off()





depths <-   read.csv("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/sampleStats_US_obtusa.csv")

dp_clean_all_merged_dp <- dp_clean_all %>%
  left_join(depths, by = c("sample" = "sampleId"))

dp_clean_all_merged_dp$standardized_depth <- dp_clean_all_merged_dp$depth/dp_clean_all_merged_dp$meanDepth

pdf("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/dp_manhattan_by_sample_100k_standardized.pdf",
    width = 25, height = 2)

for (s in unique(dp_clean_all_merged_dp$sample)) {
    sample_meanDepth <- unique(dp_clean_all_merged_dp$meanDepth[dp_clean_all_merged_dp$sample == s])
  p <- ggplot(subset(dp_clean_all_merged_dp, sample == s, meanDepth == meanDepth),
              aes(x = window_start, y = standardized_depth)) +
    geom_point(size = 0.4, alpha = 0.6, color = "darkblue") +
    facet_grid(~contig, scales = "free_x") +
    labs(
      title = paste0("Windowed heterozygosity (10kb): ", s, " | meanDepth: ", sample_meanDepth),
      x = "Genomic position (window start)",
      y = "Depth"
    ) +
    ylim(0, 2.5) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.minor = element_blank()
    )
  print(p)
}

dev.off()





dp_clean_all_merged_dp$sample <- as.factor(dp_clean_all_merged_dp$sample)
mean_profile <- dp_clean_all_merged_dp %>%
  filter(contig == "JAACYE010000002.1") %>%
  group_by(contig, window_start) %>%
  summarise(mean_std_depth = mean(standardized_depth, na.rm = TRUE),
            .groups = "drop")

pdf("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/dp_manhattan_by_sample_100k_standardized_all_chr2.pdf",
    width = 50, height = 10)

p <- ggplot(
    subset(dp_clean_all_merged_dp, contig == "JAACYE010000002.1"),
    aes(x = window_start, 
        y = standardized_depth, 
        group = sample,         # make sure each sample is its own line
        color = meanDepth)         # color by sample (clean separation)
  ) +
  geom_line() +
 geom_line(
  data = mean_profile,
  aes(x = window_start, y = mean_std_depth),
  color = "red",
  size = 1.2,
  inherit.aes = FALSE
)+
  facet_grid(~contig, scales = "free_x") +
  labs(
    title = paste("Depth profile: all"),
    x = "Genomic position (window start)",
    y = "Depth"
  ) +
  ylim(0, 2.5) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.minor = element_blank()
  )

print(p)


dev.off()





pdf("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/dp_manhattan_by_sample_100k_standardized_all_chr2_depth5.pdf",
    width = 50, height = 10)

p <- ggplot(
    subset(dp_clean_all_merged_dp, contig == "JAACYE010000002.1" & meanDepth > 5),
    aes(x = window_start, 
        y = standardized_depth, 
        group = sample,         # make sure each sample is its own line
        color = meanDepth)         # color by sample (clean separation)
  ) +
  geom_line() +
 geom_line(
  data = mean_profile,
  aes(x = window_start, y = mean_std_depth),
  color = "red",
  size = 1.2,
  inherit.aes = FALSE
)+
  facet_grid(~contig, scales = "free_x") +
  labs(
    title = paste("Depth profile: all"),
    x = "Genomic position (window start)",
    y = "Depth"
  ) +
  ylim(0, 2.5) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.minor = element_blank()
  )

print(p)


dev.off()





pdf("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/dp_manhattan_by_sample_100k_standardized_all_chr2_depthless_than_5.pdf",
    width = 50, height = 10)

p <- ggplot(
    subset(dp_clean_all_merged_dp, contig == "JAACYE010000002.1" & meanDepth < 5),
    aes(x = window_start, 
        y = standardized_depth, 
        group = sample,         # make sure each sample is its own line
        color = meanDepth)         # color by sample (clean separation)
  ) +
  geom_line() +
 geom_line(
  data = mean_profile,
  aes(x = window_start, y = mean_std_depth),
  color = "red",
  size = 1.2,
  inherit.aes = FALSE
)+
  facet_grid(~contig, scales = "free_x") +
  labs(
    title = paste("Depth profile: all"),
    x = "Genomic position (window start)",
    y = "Depth"
  ) +
  ylim(0, 2.5) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.minor = element_blank()
  )

print(p)


dev.off()





library(ggplot2)
library(dplyr)
library(readr)

# Load BED file
repeats <- read_tsv("/scratch/rjp5nc/removedups/us_dobtusa/usobtusa_repeats.bed",
                    col_names = c("chrom", "start", "end"), comment = "#")

# Add a small y-value so the repeat regions appear near the x-axis
repeats <- repeats %>%
  mutate(ymin = 0, ymax = 0.1)  # adjust ymax for visual height




png("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/repeats.png",
    width = 10000, height = 2000)
# Example genome length or windowed data for context
# Here we just plot the repeats as segments
ggplot(repeats) +
  geom_segment(aes(x = start, xend = end, y = ymin, yend = ymax),
               color = "red", size = 1) +
  facet_grid(~chrom, scales = "free_x") +
  labs(title = "Masked repeat regions",
       x = "Genomic position",
       y = "") +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
dev.off()








# Subset repeat regions for this contig
repeats_sub <- repeats %>%
  filter(chrom == "JAACYE010000002.1") %>%
  mutate(ymin = 0, ymax = 0.2)  # small track near x-axis














library(dplyr)
library(ggplot2)

# Set window size
window_size <- 100000

# Choose contig to plot
contig_to_plot <- "JAACYE010000002.1"

# Subset repeats for this contig
repeats_sub <- repeats %>%
  filter(chrom == contig_to_plot)

# Create 100kb windows along the contig
max_pos <- max(repeats_sub$end)  # or set manually if needed
windows <- data.frame(
  window_start = seq(0, max_pos, by = window_size),
  window_end = seq(window_size, max_pos + window_size, by = window_size)
)

# Count repeats per window
repeat_counts <- windows %>%
  rowwise() %>%
  mutate(num_repeats = sum(
    (repeats_sub$start < window_end) & (repeats_sub$end > window_start)
  )) %>%
  ungroup()



p <- ggplot(
    subset(dp_clean_all_merged_dp, contig == "JAACYE010000002.1" & meanDepth > 5),
    aes(x = window_start, 
        y = standardized_depth, 
        group = sample,         # make sure each sample is its own line
        color = meanDepth)         # color by sample (clean separation)
  ) +
  geom_line() +
 geom_line(
  data = mean_profile,
  aes(x = window_start, y = mean_std_depth),
  color = "red",
  size = 1.2,
  inherit.aes = FALSE
)+
  facet_grid(~contig, scales = "free_x") +
  labs(
    title = paste("Depth profile: all"),
    x = "Genomic position (window start)",
    y = "Depth"
  ) +
  ylim(0, 2.5) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.minor = element_blank()
  )


  # Plot repeat counts as a line
p2 <- ggplot(repeat_counts, aes(x = window_start, y = num_repeats)) +
  geom_line(color = "red", size = 1) +
  labs(
    title = paste("Number of repeats per 100kb window:", contig_to_plot),
    x = "Genomic position (window start)",
    y = "Number of repeats"
  ) +
  theme_bw(base_size = 12)


  combined_plot <- p / p2 + plot_layout(heights = c(3, 1))  # top plot 3x taller than bottom


repeat_prop <- windows %>%
  rowwise() %>%
  mutate(masked_bp = sum(
    pmax(
      0, 
      pmin(repeats_sub$end, window_end) - pmax(repeats_sub$start, window_start)
    )
  )) %>%
  ungroup() %>%
  mutate(prop_masked = masked_bp / window_size)

# Plot percent masked
p3 <- ggplot(repeat_prop, aes(x = window_start, y = prop_masked)) +
  geom_line(color = "darkgreen", size = 1) + ylim(0,1)+
  labs(
    title = paste("Percent of window masked by repeats:", contig_to_plot),
    x = "Genomic position (window start)",
    y = "Proportion masked"
  ) +
  theme_bw(base_size = 12)


combined_plot <- p / p2 /p3 + plot_layout(heights = c(3,1, 1))



# Save to PDF
pdf("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/dp_depth_and_repeatcounts_chr2_5andup.pdf",
    width = 60, height = 12)

print(combined_plot)

dev.off()




















library(dplyr)
library(ggplot2)

# Set window size
window_size <- 100000

# Choose contig to plot
contig_to_plot <- "JAACYE010000003.1"

# Subset repeats for this contig
repeats_sub <- repeats %>%
  filter(chrom == contig_to_plot)

# Create 100kb windows along the contig
max_pos <- max(repeats_sub$end)  # or set manually if needed
windows <- data.frame(
  window_start = seq(0, max_pos, by = window_size),
  window_end = seq(window_size, max_pos + window_size, by = window_size)
)

# Count repeats per window
repeat_counts <- windows %>%
  rowwise() %>%
  mutate(num_repeats = sum(
    (repeats_sub$start < window_end) & (repeats_sub$end > window_start)
  )) %>%
  ungroup()

mean_profile <- dp_clean_all_merged_dp %>%
  filter(contig == "JAACYE010000003.1") %>%
  group_by(contig, window_start) %>%
  summarise(mean_std_depth = mean(standardized_depth, na.rm = TRUE),
            .groups = "drop")

p4 <- ggplot(
    subset(dp_clean_all_merged_dp, contig == "JAACYE010000003.1" & meanDepth > 5),
    aes(x = window_start, 
        y = standardized_depth, 
        group = sample,         # make sure each sample is its own line
        color = meanDepth)         # color by sample (clean separation)
  ) +
  geom_line() +
 geom_line(
  data = mean_profile,
  aes(x = window_start, y = mean_std_depth),
  color = "red",
  size = 1.2,
  inherit.aes = FALSE
)+
  facet_grid(~contig, scales = "free_x") +
  labs(
    title = paste("Depth profile: all"),
    x = "Genomic position (window start)",
    y = "Depth"
  ) +
  ylim(0, 2.5) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.minor = element_blank()
  )


  # Plot repeat counts as a line
p5 <- ggplot(repeat_counts, aes(x = window_start, y = num_repeats)) +
  geom_line(color = "red", size = 1) +
  labs(
    title = paste("Number of repeats per 100kb window:", contig_to_plot),
    x = "Genomic position (window start)",
    y = "Number of repeats"
  ) +
  theme_bw(base_size = 12)


  combined_plot <- p / p2 + plot_layout(heights = c(3, 1))  # top plot 3x taller than bottom


repeat_prop <- windows %>%
  rowwise() %>%
  mutate(masked_bp = sum(
    pmax(
      0, 
      pmin(repeats_sub$end, window_end) - pmax(repeats_sub$start, window_start)
    )
  )) %>%
  ungroup() %>%
  mutate(prop_masked = masked_bp / window_size)

# Plot percent masked
p6 <- ggplot(repeat_prop, aes(x = window_start, y = prop_masked)) +
  geom_line(color = "darkgreen", size = 1) + ylim(0,1)+
  labs(
    title = paste("Percent of window masked by repeats:", contig_to_plot),
    x = "Genomic position (window start)",
    y = "Proportion masked"
  ) +
  theme_bw(base_size = 12)


combined_plot <- p4 / p5 /p6 + plot_layout(heights = c(3,1, 1))



# Save to PDF
pdf("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/dp_depth_and_repeatcounts_chr3_5andup.pdf",
    width = 60, height = 12)

print(combined_plot)

dev.off()








#combined_plot <- p + p2 + p3 + p4 + p5 + p6 +
combined_plot <- p + p4 + p2 + p5 + p3 + p6 +

  plot_layout(ncol = 2, 
              heights = c(3,1,1))  # optional: relative heights of rows

pdf("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/dp_depth_and_repeatcounts_chr3and4_5andup.pdf",
    width = 60, height = 12)

print(combined_plot)

dev.off()