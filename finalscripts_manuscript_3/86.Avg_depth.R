#ijob -A berglandlab -c10 -p standard --mem=40G

#module load gcc/11.4.0  openmpi/4.1.4 icu R/4.3.1
#R


library(ggplot2)
library(data.table)

dp_raw <- read.table("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/window_avg_depth.txt",
                     header = TRUE, sep = "\t")


dp_clean <- na.omit(dp_raw)


dp_clean$window_start <- as.numeric(dp_clean$window_start)
dp_clean$avg_depth   <- as.numeric(dp_clean$avg_depth)




dp_raw_all <- read.table("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/avg.depth.all12.txt",
                     header = FALSE, sep = "\t")

dp_clean_all <- na.omit(dp_raw_all)

colnames(dp_clean_all) <- c("contig", "window_start", "window_end", "sample", "depth")


dp_clean_all$window_start <- as.numeric(dp_clean_all$window_start)
dp_clean_all$depth   <- as.numeric(dp_clean_all$depth)


# Save plot as PNG
png("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/dp_manhattan.png",
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
png("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/dp_7_manhattan.png",
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
















pdf("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/dp_manhattan_by_sample_10000.pdf",
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



pdf("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/het_10kb_per_sample.pdf",
    width = 25, height = 2)

# Loop over samples and plot
for(s in unique(het_merged$sample)){
  p <- ggplot(subset(het_merged, sample == s),
              aes(x = win_start, y = het_prop)) +
    geom_point(size = 0.4, alpha = 0.6, color = "darkgreen") +
    facet_grid(~contig, scales = "free_x") +
    labs(
      title = paste("Sliding window heterozygosity (10kb):", s),
      x = "Genomic position (window start)",
      y = "Observed heterozygosity (Ho)"
    ) +
    ylim(0, 0.05) +  # adjust according to your het values
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.minor = element_blank()
    )
  print(p)
}

dev.off()









dp_raw_all <- read.table("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/avg.depth.all12.txt",
                     header = FALSE, sep = "\t")

dp_clean_all <- na.omit(dp_raw_all)

colnames(dp_clean_all) <- c("contig", "window_start", "window_end", "sample", "depth")


dp_clean_all$window_start <- as.numeric(dp_clean_all$window_start)
dp_clean_all$depth   <- as.numeric(dp_clean_all$depth)

pdf("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/dp_manhattan_by_sample_10000.pdf",
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






dp_raw_all <- read.table("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/avg_avgdepth.long.sorted_oneofeach_100000.txt",
                     header = FALSE, sep = " ")

dp_clean_all <- na.omit(dp_raw_all)

colnames(dp_clean_all) <- c("contig", "window_start", "window_end", "sample", "depth")


dp_clean_all$window_start <- as.numeric(dp_clean_all$window_start)
dp_clean_all$depth   <- as.numeric(dp_clean_all$depth)

pdf("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/dp_manhattan_by_sample_one_each_100000.pdf",
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








# dp_raw_all <- read.table("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/window_avg_avgdepth.long.sorted_oneofeach_100000.txt",
#                      header = TRUE, sep = "\t")

# dp_clean_all <- na.omit(dp_raw_all)
# head(dp_clean_all)


# colnames(dp_clean_all) <- c("contig", "window_start", "window_end", "depth")


# dp_clean_all$window_start <- as.numeric(dp_clean_all$window_start)
# dp_clean_all$depth   <- as.numeric(dp_clean_all$depth)



depths <-   read.csv("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/sampleStats_US_obtusa.csv")

dp_clean_all_merged_dp <- dp_clean_all %>%
  left_join(depths, by = c("sample" = "sampleId"))




dp_clean_all_merged_dp$standardized_depth <- dp_clean_all_merged_dp$depth/dp_clean_all_merged_dp$meanDepth






pdf("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/dp_manhattan_one_each_100000.pdf",
    width = 25, height = 2)

  
p <- ggplot(dp_clean_all,
              aes(x = window_start, y = depth)) +
    geom_point(size = 0.4, alpha = 0.6, color = "darkblue") +
    facet_grid(~contig, scales = "free_x") +
    labs(
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

dev.off()


















pdf("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/dp_manhattan_by_sample_10000_standardized.pdf",
    width = 25, height = 2)

for (s in unique(dp_clean_all_merged_dp$sample)) {
  p <- ggplot(subset(dp_clean_all_merged_dp, sample == s),
              aes(x = window_start, y = standardized_depth)) +
    geom_point(size = 0.4, alpha = 0.6, color = "darkblue") +
    facet_grid(~contig, scales = "free_x") +
    labs(
      title = paste("Depth profile:", s),
      x = "Genomic position (window start)",
      y = "Depth"
    ) +
    ylim(0, 10) +
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

pdf("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/dp_manhattan_by_sample_10000_standardized_all_chr2.pdf",
    width = 50, height = 50)

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
    title = paste("Depth profile: JAACYE010000002.1"),
    x = "Genomic position (window start)",
    y = "Depth"
  ) +
  ylim(0, 5) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.minor = element_blank()
  )

print(p)


dev.off()