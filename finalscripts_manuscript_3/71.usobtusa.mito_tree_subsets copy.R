##module load gcc/11.4.0  openmpi/4.1.4 icu R/4.3.1

#R

# Load required packages
library(SeqArray)
library(SNPRelate)
library(ape)
library(ggtree)
library(ggplot2)
library(dplyr)
library(tibble)
library(data.table)
library(igraph)

library(patchwork)
library(foreach)
library(lubridate)




# ---- Step 1: Open GDS file ----

metadata <- read.csv("/project/berglandlab/Robert/UKSequencing2022_2024/old_stuff/2022_2024seqmetadata20250131.csv", header = TRUE)
metadata_with_clone <- read.csv("/project/berglandlab/Robert/UKSequencing2022_2024/old_stuff/metadata_with_clone.csv", header = TRUE)

samples <- read.csv("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/Sample_ID_Species_merged_20250627.csv")

usobtusasamps <- subset(samples, Species == "Daphnia obtusa" & Continent == "NorthAmerica")
metadata_with_clone$clone <- trimws(metadata_with_clone$clone)

metadata_with_clone <- subset(metadata_with_clone, clone !="Blank")
metadata_with_clone <- subset(metadata_with_clone, clone !="BLANK")

head(metadata)
metadata$clone <- trimws(metadata$clone)
metadata <- metadata %>% 
  filter(!tolower(clone) %in% c("blank", "blanks", "na", "missing"))

seqClose(genofile)
#gds.fn <- "/scratch/rjp5nc/UK2022_2024/mito_vcf/usdobtusa_mito_genotyped.gds"             # output GDS file
gds.fn <- "/scratch/rjp5nc/UK2022_2024/mito_vcf/usdobtusa_mito_genotyped.gds"
genofile <- seqOpen(gds.fn)

seqResetFilter(genofile)

#samples_to_keep <- usobtusasamps %>% filter(Species == "Daphnia obtusa") %>% pull(Sample_ID)
samples_to_keep <- metadata_with_clone %>% filter(accuratelocation == "P58" ) %>% pull(Well)
#samples_to_keep <- metadata_with_clone %>% filter(location == "UK" | accuratelocation == "P759") %>% pull(Well)

#samples_to_keep <- metadata_with_clone$Well

unique(metadata_with_clone$accuratelocation)

seqSetFilter(genofile, sample.id = samples_to_keep)

# ---- Step 2: Filter variants with missing rate < 0.05 ----
miss_rate_per_sample <- seqMissing(genofile, per.variant = FALSE)
miss_rate_per_variant <- seqMissing(genofile, per.variant = TRUE)

sample_ids <- seqGetData(genofile, "sample.id")

valid_samples <- sample_ids[miss_rate_per_sample < 0.10]
miss_rate_per_variant <- seqMissing(genofile, per.variant=TRUE)
valid_variants <- seqGetData(genofile, "variant.id")[miss_rate_per_variant < 0.10]

final_valid_samples <- intersect(valid_samples, samples_to_keep)

seqSetFilter(genofile, sample.id = final_valid_samples)

miss_rate <- seqMissing(genofile, per.variant = TRUE)
dp <- seqGetData(genofile, "annotation/format/DP")
mean_depth <- rowMeans(dp, na.rm = TRUE)
keep <- which(miss_rate < 0.10)


# ---- Step 3: Compute IBS distance matrix ----
ibs <- snpgdsIBS(genofile, sample.id = final_valid_samples, snp.id = keep, num.thread = 4, autosome.only = FALSE)
dist_matrix <- 1 - ibs$ibs
rownames(dist_matrix) <- colnames(dist_matrix) <- ibs$sample.id
dist_matrix2 <- 1- dist_matrix

long_mat <- melt(dist_matrix2, varnames = c("CloneA", "CloneB"), value.name = "Similarity") %>%
  mutate(CloneA = as.character(CloneA),
         CloneB = as.character(CloneB))  # ensure they are character, not factor

long_mat



long_mat_plot <- ggplot(long_mat, aes(x = Similarity)) +
  geom_histogram(bins = 200) +
  theme_bw() + 
  xlab("Similarity") + ylim(0,1500) + xlim(0.5,1.05) + 
  ylab("Frequency") +
  geom_vline(xintercept = 0.975, linetype = "dashed", color = "red") 




# 2. Filter dataset
long_filt <- long_mat %>%
  filter(Similarity >= 0.95, Similarity != 1)

# 3. Build graph
g <- graph_from_data_frame(long_filt[, c("CloneA", "CloneB")], directed = FALSE)

# 4. Find connected components
comp <- components(g)

# 5. Assign group letters
group_letters <- setNames(LETTERS[comp$membership], names(comp$membership))

# 6. Add group column (both clones exist in the graph)
long_filt <- long_filt %>%
  mutate(Group = group_letters[CloneA])


head(long_filt)
unique(long_filt$Group)

long_filt_onlywell <- long_filt[,c(1,4)]

unique_clones <- long_filt_onlywell %>%
  distinct(CloneA, Group)
unique_clones
unique(long_filt$Group)

p63subset <- subset(metadata_with_clone, accuratelocation == "P58")

# ---- Step 4: Build neighbor-joining tree ----
tree <- nj(as.dist(dist_matrix))



# ---- Step 5 (Optional): Relabel tips from metadata CSV ----
# CSV should have columns: Sample_ID and Label (or Well, Clone, etc.)
tree_rooted <- root(tree, outgroup = "Rockpool1_A10", resolve.root = TRUE)

label_map <- setNames(p63subset$Well, p63subset$date)

# Apply labels (optional: fallback to Sample_ID if label is missing)
group_colors <- setNames(rainbow(length(unique(p63subset$date))), unique(p63subset$date))

# Assign color to each tip

tip_colors <- group_colors[p63subset$date]

tip_color_dt <- data.table(p63subset$date, p63subset$Well)

final_valid_samples2 <- as.data.frame(final_valid_samples)

final_valid_samples3 <- left_join(final_valid_samples2, tip_color_dt, by= c("final_valid_samples"="V2"))

tip_colors <- group_colors[final_valid_samples3$V1]


tip_color_dt <- data.table(
  sample = names(tip_colors),
  color  = as.vector(tip_colors)
)


png("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/usobtusa_mito_P58_dates.png", res = 300, width = 4000, height = 5000)
#png("/scratch/rjp5nc/UK2022_2024/mito_vcf/tree_usobtusa_circ.png", width = 1200, height = 2000)

# Plot the tree
plot.phylo(tree_rooted,
#           type = "phylogram",
           type = "fan",
           cex = 0.8,
           label.offset = 0.01,
           no.margin = TRUE,
           tip.color = tip_colors,
           main = "Rooted NJ Tree (root = P63)")

# Add legend
legend("topleft",                   # position
       legend = unique(tip_color_dt$sample),     # group names
       col = unique(tip_color_dt$color),        # matching colors
       pch = 19,                   # solid circle
       pt.cex = 1.5,               # point size
       cex = 1,                    # text size
       bty = "n",                  # no box
       title = "Sample Group")

dev.off()



