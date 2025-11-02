#ijob -A berglandlab -c2 -p standard --mem=40G
#module load gcc/11.4.0  openmpi/4.1.4 icu R/4.3.1
#R

library(tidyverse)
library(ggplot2)
library(data.table)
library(patchwork)

dp_raw <- read.table("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/all_samples_100kb_depths.tsv",
                     header = TRUE, sep = "\t")

dp_clean <- na.omit(dp_raw)


dp_clean$start <- as.numeric(dp_clean$start)
dp_clean$depth   <- as.numeric(dp_clean$depth)


depths <-   read.csv("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/sampleStats_US_obtusa.csv")

dp_clean_all_merged_dp <- dp_clean %>%
  left_join(depths, by = c("sample" = "sampleId"))

dp_clean_all_merged_dp$standardized_depth <- dp_clean_all_merged_dp$depth/dp_clean_all_merged_dp$meanDepth

dp_clean_all_merged_dp$sample <- as.factor(dp_clean_all_merged_dp$sample)

mean_profile <- dp_clean_all_merged_dp %>%
  filter(chrom == "JAACYE010000002.1") %>%
  group_by(chrom, start) %>%
  summarise(mean_std_depth = mean(standardized_depth, na.rm = TRUE),
            .groups = "drop")

p <- ggplot(
    subset(
  dp_clean_all_merged_dp,
  (chrom == "JAACYE010000002.1" & meanDepth > 5) |
  (chrom == "JAACYE010000003.1" & meanDepth > 5)),
    aes(x = start, 
        y = standardized_depth, 
        group = sample,         # make sure each sample is its own line
        color = meanDepth)         # color by sample (clean separation)
  ) +
  geom_line() +
 geom_line(
  data = mean_profile,
  aes(x = start, y = mean_std_depth),
  color = "red",
  size = 1.2,
  inherit.aes = FALSE
)+
  facet_grid(~chrom, scales = "free_x") +
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



pdf("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/dp_manhattan_100k_unfilt_bam_depth5andup.pdf",
    width = 60, height = 12)

print(p)

dev.off()



p <- ggplot(subset(
  dp_clean_all_merged_dp,
  (chrom == "JAACYE010000001.1") |
  (chrom == "JAACYE010000002.1") |
  (chrom == "JAACYE010000003.1") |
  (chrom == "JAACYE010000004.1") |
  (chrom == "JAACYE010000005.1") |
  (chrom == "JAACYE010000006.1") |
  (chrom == "JAACYE010000007.1") |
  (chrom == "JAACYE010000008.1") |
  (chrom == "JAACYE010000009.1") |
  (chrom == "JAACYE010000010.1") |
  (chrom == "JAACYE010000011.1") |
  (chrom == "JAACYE010000012.1")),
    aes(x = start, 
        y = standardized_depth, 
        group = sample,         # make sure each sample is its own line
        color = meanDepth)         # color by sample (clean separation)
  ) +
  geom_line() +
 geom_line(
  data = mean_profile,
  aes(x = start, y = mean_std_depth),
  color = "red",
  size = 1.2,
  inherit.aes = FALSE
)+
  facet_grid(~chrom, scales = "free_x") +
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



pdf("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/dp_manhattan_100k_unfilt_bam.pdf",
    width = 60, height = 12)

print(p)

dev.off()