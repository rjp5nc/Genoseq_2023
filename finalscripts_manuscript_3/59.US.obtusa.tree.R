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

# ---- Step 1: Open GDS file ----

metadata <- read.csv("/project/berglandlab/Robert/UKSequencing2022_2024/old_stuff/2022_2024seqmetadata20250131.csv", header = TRUE)
metadata_with_clone <- read.csv("/project/berglandlab/Robert/UKSequencing2022_2024/old_stuff/metadata_with_clone.csv", header = TRUE)

samples <- read.csv("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/Sample_ID_Species_merged_20250627.csv")

usobtusasamps <- subset(samples, Species == "Daphnia obtusa" & Continent == "NorthAmerica")


metadata_with_clone <- subset(metadata_with_clone, clone !="Blank" & clone !="BLANK")

head(metadata)

#gds.fn <- "/scratch/rjp5nc/UK2022_2024/mito_vcf/usdobtusa_mito_genotyped.gds"             # output GDS file
gds.fn <- "/scratch/rjp5nc/UK2022_2024/mito_vcf/usdobtusa_mito_genotyped.gds"
genofile <- seqOpen(gds.fn)

seqResetFilter(genofile)

samples_to_keep <- usobtusasamps %>% filter(Species == "Daphnia obtusa") %>% pull(Sample_ID)
#samples_to_keep <- metadata_with_clone %>% filter(accuratelocation == "P759"|accuratelocation == "Dorset") %>% pull(Well)
#samples_to_keep <- metadata_with_clone %>% filter(location == "UK" | accuratelocation == "P759") %>% pull(Well)

#samples_to_keep <- metadata_with_clone$Well

unique(metadata_with_clone$accuratelocation)

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

# 2. Filter dataset
long_filt <- long_mat %>%
  filter(Similarity >= 0.995, Similarity != 1)

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



write.csv(dist_matrix, "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/usobtusa_mito_similarity_matrix.csv")  # Save as CSV


# ---- Step 4: Build neighbor-joining tree ----
tree <- nj(as.dist(dist_matrix))

# ---- Step 5 (Optional): Relabel tips from metadata CSV ----
# CSV should have columns: Sample_ID and Label (or Well, Clone, etc.)
tree_rooted <- root(tree, outgroup = "Gilmer5_E2", resolve.root = TRUE)


label_map <- setNames(metadata_with_clone$clone, metadata_with_clone$Well)

# Apply relabeling

# Apply labels (optional: fallback to Sample_ID if label is missing)
tree$tip.label <- label_map[tree$tip.label]
tree$tip.label[is.na(tree$tip.label)] <- names(label_map)[is.na(tree$tip.label)]

group_colors <- setNames(rainbow(length(unique(metadata_with_clone$date))), unique(metadata_with_clone$date))

# Assign color to each tip

metadata_with_clone$date

tip_colors <- group_colors[metadata_with_clone$date]

tip_color_dt <- data.table(metadata_with_clone$date, metadata_with_clone$Well)

final_valid_samples2 <- as.data.frame(final_valid_samples)

final_valid_samples3 <- left_join(final_valid_samples2, tip_color_dt, by= c("final_valid_samples"="V2"))

tip_colors <- group_colors[final_valid_samples3$V1]


tip_color_dt <- data.table(
  sample = names(tip_colors),
  color  = as.vector(tip_colors)
)


png("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/usobtusa_mito_rooted_by_P759.png", res = 300, width = 4000, height = 9000)
#png("/scratch/rjp5nc/UK2022_2024/mito_vcf/tree_usobtusa_circ.png", width = 1200, height = 2000)

# Plot the tree
plot.phylo(tree_rooted,
           type = "phylogram",
#           type = "fan",
           cex = 0.8,
           label.offset = 0.01,
           no.margin = TRUE,
           tip.color = tip_colors,
           main = "Rooted NJ Tree (root = P759)")

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




png("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/usobtusa_mito_rooted_by_P759.png", res = 300, width = 4000, height = 4000)
#png("/scratch/rjp5nc/UK2022_2024/mito_vcf/tree_usobtusa_circ.png", width = 1200, height = 2000)

# Plot the tree
plot.phylo(tree_rooted,
           type = "fan",
           cex = 0.8,
           label.offset = 0.01,
           no.margin = TRUE,
           tip.color = tip_colors,
           main = "Rooted NJ Tree (root = P759)")

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





















dist_matrix <- as.matrix(dist_matrix)  # Coerce it into a matrix

library(reshape2)



dist_matrix2 <- 1 - dist_matrix



# Create a named vector for mapping
Sample<-metadata_with_clone$Well
d<-as.data.table(Sample)
#d$clone<- pca_datawithclone$clone
d$date<- metadata_with_clone$date
d$clone<- metadata_with_clone$clone

name_map <-d

name_dict <- setNames(name_map$clone, name_map$Sample)


  # Coerce it into a matrix
mode(dist_matrix2) <- "numeric"        # Ensure it's numeric

# Check for NAs and replace them (if necessary)

dist_matrix2[is.na(dist_matrix2)] <- 1  # Replace NAs with 0

# Check the dimensions to ensure it's square
if (nrow(dist_matrix2) != ncol(dist_matrix2)) {
  stop("The similarity matrix is not square!")
}



# Step 2: Replace row names in the similarity matrix
rownames(dist_matrix2) <- name_dict[rownames(dist_matrix2)]

# Step 3: Replace column names in the similarity matrix
colnames(dist_matrix2) <- name_dict[colnames(dist_matrix2)]

# subset_matrix <- sim_matrix["Elvis_03", "Elvis_03"]
# print(subset_matrix)

row_indices <- which(rownames(dist_matrix2) == "D8.4A")
col_indices <- which(colnames(dist_matrix2) == "Islands_02")

# Remove all occurrences of "B"
subset_matrix <- dist_matrix2[row_indices, col_indices]
print(subset_matrix)


pdf("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/usobtusa_mitoheatmap.pdf", width=30, height=30)


# Plot the heatmap
heatmap(dist_matrix2, 
        main = "Genotype Similarity Heatmap", 
        col = colorRampPalette(c("blue", "red"))(50), 
        scale = "none", 
        margins = c(8, 8))

dev.off()









library(phangorn)
tree <- upgma(as.dist(dist_matrix))



# ---- Step 5 (Optional): Relabel tips from metadata CSV ----
# CSV should have columns: Sample_ID and Label (or Well, Clone, etc.)


label_map <- setNames(metadata_with_clone$clone, metadata_with_clone$Well)

# Apply relabeling

# Apply labels (optional: fallback to Sample_ID if label is missing)
tree$tip.label <- label_map[tree$tip.label]
tree$tip.label[is.na(tree$tip.label)] <- names(label_map)[is.na(tree$tip.label)]




group_colors <- setNames(rainbow(length(unique(metadata_with_clone$date))), unique(metadata_with_clone$date))

# Assign color to each tip


metadata_with_clone$date

tip_colors <- group_colors[metadata_with_clone$date]

tip_color_dt <- data.table(metadata_with_clone$date, metadata_with_clone$Well)

final_valid_samples2 <- as.data.frame(final_valid_samples)

final_valid_samples3 <- left_join(final_valid_samples2, tip_color_dt, by= c("final_valid_samples"="V2"))

tip_colors <- group_colors[final_valid_samples3$V1]

tree_rooted <- root(tree, outgroup = "P759", resolve.root = TRUE)

tip_color_dt <- data.table(
  sample = names(tip_colors),
  color  = as.vector(tip_colors)
)


png("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/usobtusa_mito_all_rooted_by_P759.png", res = 300, width = 6000, height = 15000)
#png("/scratch/rjp5nc/UK2022_2024/mito_vcf/tree_usobtusa_circ.png", width = 1200, height = 2000)

# Plot the tree
plot.phylo(tree_rooted,
           type = "phylogram",
#           type = "fan",
           cex = 0.8,
           label.offset = 0.01,
           no.margin = TRUE,
           tip.color = tip_colors,
           main = "Rooted NJ Tree (root = P759)")

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




png("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/usobtusa_mitotreemitophylogram_fan_all_rooted_by_P759.png", res = 300, width = 9000, height = 9000)
#png("/scratch/rjp5nc/UK2022_2024/mito_vcf/tree_usobtusa_circ.png", width = 1200, height = 2000)

# Plot the tree
plot.phylo(tree_rooted,
           type = "fan",
           cex = 0.8,
           label.offset = 0.01,
           no.margin = TRUE,
           tip.color = tip_colors,
           main = "Rooted NJ Tree (root = P759)")

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






library(SNPRelate)


library(SeqArray)
library(SeqVarTools)

# Create SeqVarData object
svd <- SeqVarData(genofile)

# Run LD pruning using SeqVarTools
ld <- ldPruning(svd, window.size=500000, ld.threshold=0.8)

# Extract pruned SNPs
snpset <- unlist(ld)
length(snpset)
