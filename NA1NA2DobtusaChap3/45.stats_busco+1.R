#module load gcc openmpi R/4.3.1
#R
library(SeqArray)
library(SNPRelate)
library(ape)
library(ggtree)
library(ggplot2)
library(dplyr)
library(tibble)
library(data.table)

gds.fn <- "/scratch/rjp5nc/UK2022_2024/buscoanalysis/buscolifted_singlecopy.gds"
genofile <- seqOpen(gds.fn)

# Get depth matrix for all variants and samples
depth_mat <- seqGetData(genofile, "annotation/format/DP")

# For each variant, calculate proportion of samples with DP > 5
prop_high_depth <- apply(depth_mat, 2, function(dp) {
  mean(dp > 5, na.rm = TRUE)
})

# Select variants where at least 90% of samples have DP > 5
variant_pass <- prop_high_depth >= 0.99
# Convert to indices
variant_pass_idx <- which(variant_pass)

# Randomly sample 3,000 from those indices
set.seed(42)
variant_pass_idx <- sample(variant_pass_idx, 3000)

# Build logical selection from sampled indices
variant_sel <- rep(FALSE, length(prop_high_depth))
variant_sel[variant_pass_idx] <- TRUE

# Apply filter
seqSetFilter(genofile, variant.sel = variant_pass_idx, verbose = TRUE)

# Now run IBS on the subset

ibs <- snpgdsIBS(genofile, snp.id = variant_pass_idx, num.thread = 4, autosome.only = FALSE)
dist_matrix <- 1 - ibs$ibs
rownames(dist_matrix) <- colnames(dist_matrix) <- ibs$sample.id

# ---- Step 4: Build neighbor-joining tree ----
tree <- nj(as.dist(dist_matrix))





counts <- read.csv("/scratch/rjp5nc/UK2022_2024/buscoanalysis/reference_counts.csv")
species <- read.csv("/scratch/rjp5nc/UK2022_2024/buscoanalysis/species.csv")

df_summary <- counts %>%
  group_by(Filename) %>%
  mutate(Total_Count = sum(Count)) %>%  # Calculate total count per file
  mutate(Normalized_Count = Count / Total_Count)  # Divide each count by total count

df_summary <- as.data.table(df_summary)


tips_to_keep <- setdiff(
  tree$tip.label[grepl("RobertUK", tree$tip.label) | grepl("Rockpool3_A11", tree$tip.label)],
  "RobertUK_B11"
)

obtusalist <- subset(df_summary, Normalized_Count >= 0.8 & Reference == "ElvisCOI")
obtusalist$Filename <- gsub("/scratch/rjp5nc/UK2022_2024/allshortreads/counts/", "", obtusalist$Filename)
obtusalist$Filename <- gsub(".counts.txt", "", obtusalist$Filename)
obtusalist$Reference <- gsub("ElvisCOI", "D. obtusa", obtusalist$Reference)


pulexlist <- subset(df_summary, Normalized_Count >= 0.8 & Reference == "mtdna_D8_119")
pulexlist$Filename <- gsub("/scratch/rjp5nc/UK2022_2024/allshortreads/counts/", "", pulexlist$Filename)
pulexlist$Filename <- gsub(".counts.txt", "", pulexlist$Filename)
pulexlist$Reference <- gsub("mtdna_D8_119", "D. pulex", pulexlist$Reference)

ambigualist <- df_summary

ambigualist <- subset(df_summary, Normalized_Count >= 0.8 & Reference == "AF523699.1")
ambigualist$Filename <- gsub("/scratch/rjp5nc/UK2022_2024/allshortreads/counts/", "", ambigualist$Filename)
ambigualist$Filename <- gsub(".counts.txt", "", ambigualist$Filename)
ambigualist$Reference <- gsub("AF523699.1", "D. ambigua", ambigualist$Reference)

alllist <-rbind(pulexlist, obtusalist,ambigualist)
alllist$Name <- paste("mtdna.",alllist$Filename)
alllist$Name <- gsub(" ", "", alllist$Name)
alllist$Name <- factor(alllist$Name)

alllist <- alllist[,c(1,2,6)]

pca_datawithclone <- read.csv("/scratch/rjp5nc/UK2022_2024/buscoanalysis/pca_withclone.csv")

head(pca_datawithclone)

cloneidpluswell <- pca_datawithclone[,c(2,8,9,10,11,12,13)]

cloneidpluswell2 <- left_join(alllist, cloneidpluswell, by = c("Filename" = "Well"))
head(alllist)
head(cloneidpluswell)

tree_rooted <- root(tree, outgroup = "Rockpool3_A11", resolve.root = TRUE)

pruned_tree <- drop.tip(tree_rooted, setdiff(tree$tip.label, tips_to_keep))

pruned_tree_df2 <- left_join(pruned_tree, alllist, by = c("label" = "Filename"))

tree_phylo <- as.phylo(pruned_tree_df2)

# 2️⃣ Drop tips corresponding to blank clones
tips_to_drop <- cloneidpluswell2$Filename[cloneidpluswell2$clone == "BLANK"]
pruned_tree_no_blank <- drop.tip(tree_phylo, tips_to_drop)



metadata <- cloneidpluswell2 %>%
  mutate(label = Filename) %>%      # must match tree tip labels
  select(label, Reference, clone, accuratelocation)

# 4️⃣ Plot the tree with colored tip labels
buscotree <- ggtree(pruned_tree_no_blank, layout = "circular") %<+% metadata +
  geom_tiplab(aes(label = clone, color = Reference), size = 3) +
  theme_tree2()

ggsave("/scratch/rjp5nc/UK2022_2024/buscoanalysis/busco_tree_UK_3000.pdf", plot = buscotree, width = 20, height = 20, dpi = 300, limitsize = FALSE)


buscotree <- ggtree(pruned_tree_no_blank, layout = "circular") %<+% metadata +
  geom_tiplab(aes(label = clone, color = accuratelocation), size = 3) +
  theme_tree2()

ggsave("/scratch/rjp5nc/UK2022_2024/buscoanalysis/busco_tree_UK_3000_accuratelocation.pdf", plot = buscotree, width = 20, height = 20, dpi = 300, limitsize = FALSE)






allmetadata <- read.csv("/scratch/rjp5nc/All_clones_metadata.csv")

# Select all North American D. obtusa tips
north_america_tips <- allmetadata$Sample_ID[
  allmetadata$Species == "Daphnia obtusa" & allmetadata$Continent == "NorthAmerica"
]

# Combine with the previous selection, excluding RobertUK_B11
tips_to_keep <- setdiff(
  c(
    tree$tip.label[grepl("Rockpool3_A11", tree$tip.label)],
    north_america_tips
  ),
  "RobertUK_B11")

tree_rooted <- root(tree, outgroup = "Rockpool3_A11", resolve.root = TRUE)

pruned_tree <- drop.tip(tree_rooted, setdiff(tree$tip.label, tips_to_keep))

tips_to_drop <- c(
  cloneidpluswell2$Filename[cloneidpluswell2$clone == "BLANK"],
  cloneidpluswell2$Filename[cloneidpluswell2$clone == "Blank"]
)
pruned_tree_no_blank <- drop.tip(pruned_tree, tips_to_drop)


allmetadata2 <- cloneidpluswell %>%
  mutate(label = Well) %>%      # must match tree tip labels
  select(Well, clone, date, accuratelocation)

buscotree <- ggtree(pruned_tree_no_blank, layout = "circular") %<+% allmetadata2 +
  geom_tiplab(aes(label = accuratelocation, color = date), size = 3) +
  theme_tree2()

ggsave("/scratch/rjp5nc/UK2022_2024/buscoanalysis/busco_tree_NA_3000.pdf", plot = buscotree, width = 20, height = 20, dpi = 300, limitsize = FALSE)
















allmetadata <- read.csv("/scratch/rjp5nc/All_clones_metadata.csv")

# Select all North American D. obtusa tips
obtusa_tips <- allmetadata$Sample_ID[
  allmetadata$Species == "Daphnia obtusa" & allmetadata$Continent == "NorthAmerica" | allmetadata$Species == "Daphnia obtusa" & allmetadata$Continent == "Europe"
]

# Combine with the previous selection, excluding RobertUK_B11
tips_to_keep <- setdiff(
  c(
    tree$tip.label[grepl("Rockpool3_A11", tree$tip.label)],
    obtusa_tips
  ),
  "RobertUK_B11")

tree_rooted <- root(tree, outgroup = "Rockpool3_A11", resolve.root = TRUE)

pruned_tree <- drop.tip(tree_rooted, setdiff(tree$tip.label, tips_to_keep))

tips_to_drop <- c(
  cloneidpluswell2$Filename[cloneidpluswell2$clone == "BLANK"],
  tree$tip.label[grepl("March20_2018_DBunk_18", tree$tip.label)], 
  tree$tip.label[grepl("March20_2018_DBunk_38", tree$tip.label)],
  tree$tip.label[grepl("March20_2018_DBunk_10", tree$tip.label)],
  cloneidpluswell2$Filename[cloneidpluswell2$clone == "Blank"]
)

pruned_tree_no_blank <- drop.tip(pruned_tree, tips_to_drop)


allmetadata2 <- cloneidpluswell %>%
  mutate(label = Well) %>%      # must match tree tip labels
  select(Well, clone, date, accuratelocation)

buscotree <- ggtree(pruned_tree_no_blank, layout = "circular") %<+% allmetadata +
  geom_tiplab(aes(label = label, color = Continent), size = 3) +
  theme_tree2()

ggsave("/scratch/rjp5nc/UK2022_2024/buscoanalysis/busco_tree_obtusa_3000.pdf", plot = buscotree, width = 20, height = 20, dpi = 300, limitsize = FALSE)






#snprelate for K0 by IBS