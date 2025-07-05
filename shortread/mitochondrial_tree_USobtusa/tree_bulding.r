
#module load gcc/11.4.0  openmpi/4.1.4 icu R/4.3.1

library(vcfR)
library(adegenet)
library(ape)
library(dplyr)
library(proxy)
library(ggtree)

vcf <- read.vcfR("/scratch/rjp5nc/UK2022_2024/mito_vcf/mito_lifted_all_missingasref.vcf.gz")

# Step 2: Extract Genotype Data
genotypes <- extract.gt(vcf)

# Step 3: Convert Genotype Data to Genind Object
genind_obj <- vcfR2genind(vcf)


metadata <- read.csv("/project/berglandlab/Robert/UKSequencing2022_2024/old_stuff/2022_2024seqmetadata20250131.csv", header = TRUE)
metadata_with_clone <- read.csv("/project/berglandlab/Robert/UKSequencing2022_2024/old_stuff/metadata_with_clone.csv", header = TRUE)

head(metadata)

samples_to_keep <- metadata_with_clone %>% filter(accuratelocation == "P63"|accuratelocation == "P62"|accuratelocation == "PBO66"|accuratelocation == "P58"|accuratelocation == "Gilmer"|accuratelocation == "P759") %>% pull(Well)
  
#samples_to_keep <- metadata_with_clone %>% pull(Sample)
#
#|accuratelocation == "P63"|accuratelocation == "P62"|accuratelocation == "PBO66"|accuratelocation == "P58"|clone=="simo"
# & Plate != "simo"

## Step 6: Subset Genind Object to Keep Only Selected Samples
genind_obj2 <- genind_obj[indNames(genind_obj) %in% samples_to_keep, ]

gen_matrix <- tab(genind_obj2)

# Step 4: Compute Distance Matrix

dist_matrix <- dist(gen_matrix)


# Step 2: Make sure the Sample_ID matches the tree tip labels
# Check this with: `treemito$tip.label`
# Then relabel:

treemito <- njs(dist_matrix)

# Clean whitespaces
treemito$tip.label <- trimws(treemito$tip.label)
metadata_with_clone$Well <- trimws(metadata_with_clone$Well)
metadata_with_clone$clone <- trimws(metadata_with_clone$clone)

# Make the mapping from Well to clone
label_map <- setNames(metadata_with_clone$clone, metadata_with_clone$Well)

# Apply relabeling
treemito$tip.label <- label_map[treemito$tip.label]

# Confirm relabeling worked
head(treemito$tip.label)
table(is.na(treemito$tip.label))  


png("/scratch/rjp5nc/UK2022_2024/mito_vcf/treemito_circular_labeled.png", width = 1200, height = 1200)
plot.phylo(treemito,
           type = "fan",           # circular layout
           cex = 0.8,              # tip label size
           label.offset = 0.01,    # spacing for labels
           no.margin = TRUE,
           main = "Circular Neighbor-Joining Tree (mito)")
dev.off()










#sim_matrix <- 1 / (1 + (dist_matrix))  # Convert to similarity (scaled 0-1)

#genind_obj3<- tab(genind_obj2)
#sim_matrix <- as.matrix(simil(genind_obj, method = "Jaccard"))

heatmap(sim_matrix)  # Quick visualization




genind_data <- tab(genind_obj2)  # Convert genind object to genotype table
sim_matrix <- simil(genind_data, method = "Jaccard")

# Optionally, convert to a numeric matrix for further processing
sim_matrix <- as.matrix(sim_matrix)

write.csv(sim_matrix, "/scratch/rjp5nc/UK2022_2024/mito_vcf/fullmito_similarity_matrixobtusa.csv")  # Save as CSV

















# Load required packages
library(SeqArray)
library(SNPRelate)
library(ape)
library(ggtree)
library(ggplot2)
library(dplyr)
library(tibble)

# ---- Step 1: Open GDS file ----



metadata <- read.csv("/project/berglandlab/Robert/UKSequencing2022_2024/old_stuff/2022_2024seqmetadata20250131.csv", header = TRUE)
metadata_with_clone <- read.csv("/project/berglandlab/Robert/UKSequencing2022_2024/old_stuff/metadata_with_clone.csv", header = TRUE)


metadata_with_clone <- subset(metadata_with_clone, clone !="Blank" & clone !="BLANK")

head(metadata)

seqResetFilter(genofile)


samples_to_keep <- metadata_with_clone %>% filter(accuratelocation == "P63"|accuratelocation == "P62"|accuratelocation == "PBO66"|accuratelocation == "P58"|accuratelocation == "Gilmer"|accuratelocation == "P759"|accuratelocation == "Dorset") %>% pull(Well)
samples_to_keep <- metadata_with_clone %>% filter(accuratelocation == "P759"|accuratelocation == "Dorset") %>% pull(Well)

samples_to_keep <- metadata_with_clone$Well

unique(metadata_with_clone$accuratelocation)

# ---- Step 2: Filter variants with missing rate < 0.05 ----
miss_rate_per_sample <- seqMissing(genofile, per.variant = FALSE)
miss_rate_per_variant <- seqMissing(genofile, per.variant = TRUE)
sample_ids <- seqGetData(genofile, "sample.id")
valid_samples <- sample_ids[miss_rate_per_sample < 0.05]
miss_rate_per_variant <- seqMissing(genofile, per.variant=TRUE)
valid_variants <- seqGetData(genofile, "variant.id")[miss_rate_per_variant < 0.05]




final_valid_samples <- intersect(valid_samples, samples_to_keep)


seqSetFilter(genofile, sample.id = final_valid_samples)

miss_rate <- seqMissing(genofile, per.variant = TRUE)
dp <- seqGetData(genofile, "annotation/format/DP")
mean_depth <- rowMeans(dp, na.rm = TRUE)
keep <- which(miss_rate < 0.05 & mean_depth > 10)


# ---- Step 3: Compute IBS distance matrix ----
ibs <- snpgdsIBS(genofile, sample.id = final_valid_samples, snp.id = keep, num.thread = 4, autosome.only = FALSE)
dist_matrix <- 1 - ibs$ibs
rownames(dist_matrix) <- colnames(dist_matrix) <- ibs$sample.id

# ---- Step 4: Build neighbor-joining tree ----
tree <- nj(as.dist(dist_matrix))

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


png("/scratch/rjp5nc/UK2022_2024/mito_vcf/treemitophylogram_circular_rooted_by_P759.png", width = 1200, height = 2000)

# Plot the tree
plot.phylo(tree_rooted,
           type = "phylogram",
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




