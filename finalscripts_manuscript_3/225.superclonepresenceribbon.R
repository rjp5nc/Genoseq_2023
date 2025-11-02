#ijob -A berglandlab -c10 -p standard --mem=40G
#module load gcc/11.4.0  openmpi/4.1.4 icu R/4.3.1
#R


library(ggplot2)
library(ggalluvial)
library(dplyr)
library(viridis)
library(dplyr)
library(data.table)

metadata_with_clone <- read.csv("/project/berglandlab/Robert/UKSequencing2022_2024/old_stuff/metadata_with_clone.csv", header = TRUE)
metadata_with_clone <- subset(metadata_with_clone, clone !="Blank" & clone !="BLANK")
mitotypes <- read.csv("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/mito_types.csv")
genomic_types <- read.csv("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/genomic_types.csv")



pops <- read.table("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/pops_fixed.txt",
                   header = FALSE, stringsAsFactors = FALSE)

genomic_ids <- genomic_types$CloneA
pop_ids <- pops$V1

# find those in pops but not in genomic_types
pops_not_in_genomic <- setdiff(pop_ids, genomic_ids)

# make a dataframe for these, with Group = "OO"
new_rows <- data.frame(
  CloneA = pops_not_in_genomic,
  Group = paste0("OO_", seq_along(pops_not_in_genomic)),
  stringsAsFactors = FALSE
)

# bind to genomic_types
genomic_types_extended <- rbind(
  genomic_types[, c("CloneA", "Group")],  # keep only relevant cols
  new_rows
)

genomic_types2 <- genomic_types_extended %>%
  left_join(metadata_with_clone, by = c("CloneA" = "Well")) %>%
  left_join(mitotypes, by = c("CloneA" = "CloneA"))

genomic_types2 <- subset(genomic_types2, date != "2023" & date != "2024")


genomic_types_filtered <- genomic_types2 %>%
  select(CloneA, date, Group.x, accuratelocation, Group.y) %>%
  rename(superclone = Group.x,
         mitotype = Group.y)


genomic_types_clean <- genomic_types_filtered %>%
  na.omit()

# Count (sum) by superclone and date
superclone_counts <- genomic_types_clean %>%
  group_by(superclone, accuratelocation, date) %>%
  summarise(count = n(), .groups = "drop")

# Count (sum) by mitotype and date
mitotype_counts <- genomic_types_clean %>%
  group_by(mitotype, accuratelocation, date) %>%
  summarise(count = n(), .groups = "drop")


mitotype_counts <- as.data.frame(mitotype_counts)

mitotype_counts$date <- as.Date(mitotype_counts$date, format = "%m/%d/%Y")
mitotype_counts$date[is.na(mitotype_counts$date)] <- as.Date(paste0(mitotype_counts$date[is.na(mitotype_counts$date)], "-01-01"))

# Reorder factor levels by date
mitotype_counts$date <- factor(mitotype_counts$date, levels = unique(mitotype_counts$date[order(mitotype_counts$date)]))



png("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/mitotype_counts.png", res = 600,
    width = 10000, height = 5000)

# Plot
ggplot(mitotype_counts,
       aes(x = date, 
           stratum = mitotype, 
           alluvium = mitotype, 
           y = count,
           fill = mitotype,
           label = mitotype)) +
  geom_flow(alpha = 0.5, width = 0.3) +
  geom_stratum(width = 0.3, color = "grey30") +
  geom_text(stat = "stratum", size = 3) +
  scale_fill_viridis_d(option = "C") +
  theme_bw() + facet_grid(~accuratelocation, scales = "free_x")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "mitotype Transitions Over Time",
       y = "Count",
       x = "Timepoint")

dev.off()









superclone_counts <- as.data.frame(superclone_counts)

superclone_counts$date <- as.Date(superclone_counts$date, format = "%m/%d/%Y")
superclone_counts$date[is.na(superclone_counts$date)] <- as.Date(paste0(superclone_counts$date[is.na(superclone_counts$date)], "-01-01"))

# Reorder factor levels by date
superclone_counts$date <- factor(superclone_counts$date, levels = unique(superclone_counts$date[order(superclone_counts$date)]))



# Identify OO* superclones
oo_clones <- unique(superclone_counts$superclone[grepl("^OO", superclone_counts$superclone)])
other_clones <- setdiff(unique(superclone_counts$superclone), oo_clones)

# Assign colors: gray for OO*, viridis for others
fill_colors <- c(
  setNames(rep("grey80", length(oo_clones)), oo_clones),
  setNames(viridis(length(other_clones), option = "C"), other_clones)
)

png("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/superclone_counts.png", res = 600,
    width = 10000, height = 5000)

ggplot(superclone_counts,
       aes(x = date,
           stratum = superclone,
           alluvium = superclone,
           y = count,
           fill = superclone,
           label = superclone)) +
  geom_flow(alpha = 0.5, width = 0.3) +
  geom_stratum(width = 0.3, color = "grey30") +
  geom_text(stat = "stratum", size = 3) +
  scale_fill_manual(
    values = fill_colors,
    breaks = other_clones  # only show non-OO* in legend
  ) +
  theme_bw() +
  facet_grid(~accuratelocation, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Superclones Over Time",
       y = "Count",
       x = "Timepoint")

dev.off()























# Count (sum) by superclone and date
superclone_counts <- genomic_types_clean %>%
  group_by(superclone, accuratelocation) %>%
  summarise(count = n(), .groups = "drop")

superclone_counts <- as.data.frame(superclone_counts)

superclone_counts$date <- as.Date(superclone_counts$date, format = "%m/%d/%Y")
superclone_counts$date[is.na(superclone_counts$date)] <- as.Date(paste0(superclone_counts$date[is.na(superclone_counts$date)], "-01-01"))

# Reorder factor levels by date
superclone_counts$date <- factor(superclone_counts$date, levels = unique(superclone_counts$date[order(superclone_counts$date)]))



# Identify OO* superclones
oo_clones <- unique(superclone_counts$superclone[grepl("^OO", superclone_counts$superclone)])
other_clones <- setdiff(unique(superclone_counts$superclone), oo_clones)

# Assign colors: gray for OO*, viridis for others
fill_colors <- c(
  setNames(rep("grey80", length(oo_clones)), oo_clones),
  setNames(viridis(length(other_clones), option = "C"), other_clones)
)

png("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/superclone_counts_overponds.png",
    width = 1000, height = 500)

ggplot(superclone_counts,
       aes(x = accuratelocation,
           stratum = superclone,
           alluvium = superclone,
           y = count,
           fill = superclone,
           label = superclone)) +
  geom_flow(alpha = 0.5, width = 0.3) +
  geom_stratum(width = 0.3, color = "grey30") +
  geom_text(stat = "stratum", size = 3) +
  scale_fill_manual(
    values = fill_colors,
    breaks = other_clones  # only show non-OO* in legend
  ) +
  theme_bw() +
  facet_grid(~accuratelocation, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Superclones Over Time",
       y = "Count",
       x = "Timepoint")

dev.off()









png("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/superclone_counts_overponds.png", res = 600,
    width = 5000, height = 5000)

ggplot(superclone_counts,
       aes(x = accuratelocation,
           stratum = superclone,
           alluvium = superclone,
           y = count,
           fill = superclone,
           label = superclone)) +
  geom_flow(alpha = 0.5, width = 0.3) +
  geom_stratum(width = 0.3, color = "grey30") +
  geom_text(stat = "stratum", size = 3) +
  scale_fill_manual(
    values = fill_colors,
    breaks = other_clones  # only show non-OO* in legend
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Superclones Over Time",
       y = "Count",
       x = "Timepoint")

dev.off()



















# bcftools query -f '[%SAMPLE\t%GT\n]' /scratch/rjp5nc/UK2022_2024/mito_vcf/merged_gvcf/usdobtusa_mito_combined.g.vcf.gz | \
#  awk '{print $1,$2}' > /scratch/rjp5nc/UK2022_2024/mito_vcf/merged_gvcf/usdobtusa_mito_combined.genotypes.tsv


library(vcfR)

# Read your VCF file
vcf <- read.vcfR("/scratch/rjp5nc/snapp5/usdobtusa_mito_biallelic.vcf.gz")
samples <- colnames(vcf@gt)[-1]  # first column is "FORMAT", remove it
samples
#14601 sites

samplesdf <- as.data.frame(samples)


samplesdf2 <- samplesdf %>%
  left_join(genomic_types_clean, by = c("samples" = "CloneA"))


ind1 <- "Gilmer5_F7"
ind2 <- "Gilmer5_G6"

# Function to convert genotypes to numeric
geno_numeric <- function(gt) {
  gt[gt %in% c("./.", ".")] <- NA
  num <- rep(NA, length(gt))
  num[gt %in% c("0/0", "0|0")] <- 0
  num[gt %in% c("0/1", "1/0", "0|1", "1|0")] <- 1
  num[gt %in% c("1/1", "1|1")] <- 2
  return(num)
}


g1 <- gt[, ind1]
g2 <- gt[, ind2]

# Convert genotypes
g1_num <- geno_numeric(g1)
g2_num <- geno_numeric(g2)

# Keep only positions present in both individuals
valid_sites <- !is.na(g1_num) & !is.na(g2_num)

# Compute dxy
dxy <- mean(abs(g1_num[valid_sites] - g2_num[valid_sites]))
dxy
mu <- 1.75e-7          # mutation rate per site per generation
gens_per_year <- 6     # adjust to your species if known

genome_size <- 14601   # or total callable sites in your genome
dxy_genome <- dxy * (length(g1_num) / genome_size)
years_divergence <- dxy_genome / (2 * mu) * gens_per_year
years_divergence