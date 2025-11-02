#ijob -A berglandlab -c10 -p standard --mem=100G
#module load gcc/11.4.0  openmpi/4.1.4 icu R/4.3.1
#R

library(data.table)
pi_raw <- read.table("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/pixy_allusobtusa/pixy_10000__pi.txt",header = TRUE,sep = "\t")

library(ggplot2)

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
png("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/pi_manhattan_onlyobtusa.png",
    width = 10000, height = 5000, res = 300)

ggplot(pi_clean, aes(x = window_pos_1, y = avg_pi)) +
  geom_point(size = 0.5, alpha = 0.6, color = "darkblue") +
  facet_grid(pop ~ chromosome, scales = "free_x", space = "free_x") +
  labs(
    title = "Manhattan-style pi Plot",
    x = "Genomic position (window start)",
    y = "pi"
  ) + ylim(0,0.05)+
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
png("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/pi_manhattan_onlyobtusa_line.png",
    width = 10000, height = 5000, res = 300)

ggplot(pi_clean, aes(x = window_pos_1, y = avg_pi)) +
  geom_line() +
  facet_grid(pop ~ chromosome, scales = "free_x", space = "free_x") +
  labs(
    title = "Manhattan-style pi Plot",
    x = "Genomic position (window start)",
    y = "pi"
  ) + ylim(0,0.05)+
  theme_bw(base_size = 12) +
  theme(
    strip.text.x = element_text(size = 10, angle = 90),
    strip.text.y = element_text(size = 8),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.minor = element_blank()
  )

dev.off()




png("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/pi_manhattan_onlyobtusa_boxplot_chrom.png",
    width = 5000, height = 2500, res = 300)

ggplot(pi_clean, aes(x = pop, y = avg_pi)) +
  geom_boxplot() +   ylim(0,0.05)+ facet_grid(~chromosome)+
  theme_bw() 

dev.off()



png("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/pi_manhattan_onlyobtusa_boxplot.png",
    width = 5000, height = 2500, res = 300)

ggplot(pi_clean, aes(x = pop, y = avg_pi)) +
  geom_boxplot() +   ylim(0,0.05)+ 
  theme_bw() 

dev.off()