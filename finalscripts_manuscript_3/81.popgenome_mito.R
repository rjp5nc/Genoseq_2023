#ijob -A berglandlab -c10 -p standard --mem=40G

#module load gcc/11.4.0  openmpi/4.1.4 icu R/4.3.1

#R

#install.packages(c("usethis","pkgdown","gh","rcmdcheck"), repos="https://cloud.r-project.org")


#install.packages("devtools")
#library(devtools)

#devtools::install_github("pievos101/PopGenome")


library(PopGenome)

# 1. Load the VCF
vcf_file <- "/scratch/rjp5nc/UK2022_2024/mito_vcf/usdobtusa_mito_genotyped_subset.vcf.gz"
genome <- readVCF(
  vcf_file,
  numcols = 10000,
  tid = "CM028013.1",   # <-- replace with your contig/chromosome name from the VCF header
  frompos = 1,
  topos = 999999999,   # just make this very large so it covers the whole contig
  include.unknown = TRUE
)


pops <- read.table("/scratch/rjp5nc/UK2022_2024/mito_vcf/pops_fixed_filtered.txt",
                   header = FALSE, stringsAsFactors = FALSE)
colnames(pops) <- c("sample", "population")

to_remove <- pops$population == "PBO66"

# Keep only non-PBO66 samples and their population labels
pops$sample <- pops$sample[!to_remove]
pops$population <- pops$population[!to_remove]

# Split into list of populations
pop_list <- split(pops$sample, pops$population)

# Set populations
genome <- set.populations(genome, pop_list)

# Set Gilmer as outgroup
genome <- set.outgroup(genome, pop_list$Gilmer, diploid = FALSE)

# Calculate π (within-pop diversity)
genome <- diversity.stats(genome)
pi_vals <- genome@nuc.diversity.within

# Calculate Dxy (between populations)
genome <- diversity.stats.between(genome)
dxy_matrix <- genome@nuc.diversity.between

# Show results
pi_vals
dxy_matrix

colnames(pi_vals) <- names(pop_list)

# Rename Dxy results
colnames(dxy_matrix) <- combn(names(pop_list), 2, FUN = paste, collapse = "/")


library(reshape2)
library(ggplot2)

# Convert your Dxy row into a full square matrix
# (PopGenome outputs in long "pairwise" style)
popslist <- names(pop_list)

# Create empty square matrix
dxy_full <- matrix(NA, nrow = length(popslist), ncol = length(popslist),
                   dimnames = list(popslist, popslist))

# Fill in the pairwise values
pairs <- combn(popslist, 2, simplify = FALSE)
vals <- as.numeric(dxy_matrix[1,])  # Dxy values
for (i in seq_along(pairs)) {
  a <- pairs[[i]][1]
  b <- pairs[[i]][2]
  dxy_full[a, b] <- vals[i]
  dxy_full[b, a] <- vals[i]  # mirror
}
diag(dxy_full) <- 0  # self comparisons

# Convert to long format for ggplot
dxy_long <- melt(dxy_full, varnames = c("Pop1", "Pop2"), value.name = "Dxy")


png("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/dxy.png",
    width = 2000, height = 2000, res = 300)

# Plot heatmap
ggplot(dxy_long, aes(x = Pop1, y = Pop2, fill = Dxy)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(option = "plasma", na.value = "white") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Pairwise Dxy between Populations",
       fill = "Dxy")

dev.off()
