#ijob -A berglandlab -c2 -p standard --mem=40G
#module load gcc/11.4.0  openmpi/4.1.4 icu R/4.3.1
#R


library(tidyverse)
library(gridExtra)
library(ggplot2)
# Read the vcftools --het output



het_data <- read.table("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/eudobtusa_indv/euobtusa_het_100kb.txt", header = FALSE, stringsAsFactors = FALSE)

colnames(het_data) <- c("contig", "win_start", "win_end", "sample", "het_prop")

# List all files


dp_raw_all <- read.table("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/eudobtusa_indv/avgdepth_all_unfilt_100k.txt",
                     header = FALSE, sep = " ")

dp_clean_all <- na.omit(dp_raw_all)

colnames(dp_clean_all) <- c("contig", "win_start", "win_end", "sample", "depth")


dp_clean_all$win_start <- as.numeric(dp_clean_all$win_start)
dp_clean_all$depth   <- as.numeric(dp_clean_all$depth)




het_merged <- het_data %>%
  left_join(dp_clean_all, by = c("contig", "win_start", "win_end", "sample"))

het_merged <- subset(het_merged, depth >= 6)



# %>%
#   left_join(genomic_types, by = c("sample" = "CloneA")) %>%
#   left_join(depths, by = c("sample" = "sampleId"))

# Plot
pdf("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/eudobtusa_indv/het_window_all_100000.pdf",
    width = 25, height = 2)

for(s in unique(het_merged$sample)){
  sample_info <- het_merged %>% filter(sample == s) %>% slice(1)

  p <- ggplot(subset(het_merged, sample == s),
              aes(x = win_start, y = het_prop)) +
    geom_point(size = 0.4, alpha = 0.6, color = "darkgreen") +
    facet_grid(~contig, scales = "free_x", space="free_x") +
    labs(
      title = paste0("Windowed heterozygosity (100kb): ", s),
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





