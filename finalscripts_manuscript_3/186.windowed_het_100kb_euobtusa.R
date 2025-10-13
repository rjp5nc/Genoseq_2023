#ijob -A berglandlab -c10 -p standard --mem=40G
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





# Plot
pdf("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/het_sliding_window_Gilmer_100000.pdf",
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
      title = paste0("Sliding window heterozygosity (10kb): ", s,
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