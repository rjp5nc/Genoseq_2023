#ijob -A berglandlab -c10 -p standard --mem=40G
#module load gcc/11.4.0  openmpi/4.1.4 icu R/4.3.1
#R


library(tidyverse)
library(gridExtra)
library(ggplot2)
# Read the vcftools --het output
het <- read.table("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/heterozygote/sample_het.het", header = TRUE)
genomic_types <- read.csv("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/genomic_types.csv")


het <- het %>%
  mutate(Ho = 1 - (O.HOM. / N_SITES))

# Merge by sample name
# Note: het$INDV = CloneA column in genomic_types
het_merged <- het %>%
  left_join(genomic_types, by = c("INDV" = "CloneA"))

# Preview merged data
head(het_merged)

# Order samples by Group
het_merged <- het_merged %>%
  arrange(Group, INDV) %>%
  mutate(INDV = factor(INDV, levels = INDV))  # preserve order

# Identify positions for group labels (middle of each group)
group_positions <- het_merged %>%
  group_by(Group) %>%
  summarise(x = mean(as.numeric(INDV)), .groups = "drop")

# Plot bars
p1 <- ggplot(het_merged, aes(x = INDV, y = Ho, fill = Group)) +
  geom_bar(stat = "identity") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_blank(),   # hide default x-axis labels
    axis.ticks.x = element_blank()
  ) +
  labs(
    title = "Per-sample Observed Heterozygosity",
    x = "Group",
    y = "Observed heterozygosity (Ho)",
    fill = "Group"
  ) +
  # Add custom group labels at the bottom
  geom_text(data = group_positions,
            aes(x = x, y = -0.02, label = Group),
            inherit.aes = FALSE, angle = 0, vjust = 1, size = 5)




# # Plot per-sample heterozygosity with x-axis sorted by Group
# p1 <- ggplot(het_merged, aes(x = INDV, y = Ho, fill = Group)) +
#   geom_bar(stat = "identity") +
#   theme_minimal(base_size = 14) +
#   theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
#   labs(
#     title = "Per-sample Observed Heterozygosity",
#     x = "Sample",
#     y = "Observed heterozygosity (Ho)",
#     fill = "Group"
#   )


# Compute per-group average heterozygosity
het_group <- het_merged %>%
  group_by(Group) %>%
  summarise(mean_Ho = mean(Ho, na.rm = TRUE))

# Plot 2: per-group heterozygosity
p2 <- ggplot(het_group, aes(x = Group, y = mean_Ho, fill = Group)) +
  geom_bar(stat = "identity") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Mean Observed Heterozygosity per Group",
    x = "Group",
    y = "Mean Ho",
    fill = "Group"
  )

# Save plots
png("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/heterozygote/heterozygote_sites.png",
    width = 800, height = 1000)

grid.arrange(p1, p2, ncol = 1)

dev.off()






# Merge heterozygosity and depth info
het_depth <- het_merged %>%
  left_join(depths, by = c("INDV" = "sampleId"))

# Order samples by Group or meanDepth
het_depth <- het_depth %>%
  arrange(Group, INDV) %>%
  mutate(INDV = factor(INDV, levels = INDV))  # preserve order

# Optional: compute group midpoints for custom x-axis labels
group_positions <- het_depth %>%
  group_by(Group) %>%
  summarise(x = mean(as.numeric(INDV)), .groups = "drop")

# Plot bars colored by meanDepth
p_depth <- ggplot(het_depth, aes(x = INDV, y = Ho, fill = meanDepth)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis_c(option = "viridis") +   # nice continuous color scale
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_blank(),  # hide individual sample names
    axis.ticks.x = element_blank()
  ) +
  labs(
    title = "Per-sample Observed Heterozygosity Colored by meanDepth",
    x = "Group",
    y = "Observed heterozygosity (Ho)",
    fill = "meanDepth"
  ) +
  geom_text(data = group_positions,
            aes(x = x, y = -0.02, label = Group),
            inherit.aes = FALSE, angle = 0, vjust = 1, size = 5)

# Save to PNG
png("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/heterozygote/heterozygosity_by_depth.png",
    width = 800, height = 600)
print(p_depth)
dev.off()




het_merged <- het_merged %>%
  arrange(Group, F) %>%
  mutate(INDV = factor(INDV, levels = INDV))

x_labels <- het_merged %>%
  arrange(Group) %>%                # keep bars sorted by Group
  mutate(label = ifelse(!duplicated(Group), as.character(Group), "")) %>%
  pull(label)

# Plot
png("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/heterozygote/F_per_sample.png",
    width = 800, height = 600)

ggplot(het_merged %>% arrange(Group), aes(x = INDV, y = F, fill = Group)) +
  geom_bar(stat = "identity") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_x_discrete(labels = x_labels) +
  labs(
    title = "Per-sample Inbreeding Coefficient (F)",
    x = "Sample",
    y = "F (inbreeding coefficient)",
    fill = "Group"
  )

dev.off()


# Merge heterozygosity and depth info
het_depths <- het_merged %>%
  left_join(depths, by = c("INDV" = "sampleId"))

# Arrange by Group
het_depths <- het_depths %>%
  arrange(Group, F) %>%
  mutate(xpos = 1:n())

# Determine x-axis labels: first sample of each group
group_labels <- het_depths %>%
  group_by(Group) %>%
  slice(1) %>%      # take the first sample per group
  ungroup() %>%
  select(xpos, Group)

# Plot
png("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/heterozygote/F_per_sample_depth.png",
    width = 1000, height = 600)

ggplot(het_depths, aes(x = xpos, y = F, fill = meanDepth)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis_c(option = "C") +
  scale_x_continuous(
    breaks = group_labels$xpos,
    labels = group_labels$Group
  ) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(
    title = "Per-sample Inbreeding Coefficient (F) colored by meanDepth",
    x = "Group",
    y = "F (inbreeding coefficient)",
    fill = "meanDepth"
  )

dev.off()