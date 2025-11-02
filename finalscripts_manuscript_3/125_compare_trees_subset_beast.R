##module load gcc/11.4.0  openmpi/4.1.4 icu R/4.3.1

#R

# Load libraries
library(ape)
library(phytools)
library(dendextend)

# Load required packages
library(SeqArray)
library(SNPRelate)
library(ggtree)
library(ggplot2)
library(dplyr)
library(tibble)
library(data.table)
library(igraph)
library(SeqVarTools)
library(tidyverse)

library(patchwork)
library(foreach)
library(lubridate)

# -----------------------------
# 1. Load NJ trees
# -----------------------------

tree1_file <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/usdobtusa_tree_genomic_only2.nwk"
tree2 <- read.nexus("/scratch/rjp5nc/snapp5/snapp2.tree")

tree1 <- read.tree(tree1_file)
depths <-   read.csv("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/sampleStats_US_obtusa.csv")
genomic_types <- read.csv("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/genomic_types.csv")



genomic_types2 <- genomic_types %>%
  left_join(depths, by = c("CloneA" = "sampleId"))



genomic_types3 <- genomic_types2 %>%
  filter(CloneA %in% tree1$tip.label)

# For each Group, keep the CloneA with highest meanDepth
best_samples <- genomic_types3 %>%
  group_by(Group) %>%
  slice_max(meanDepth, n = 1, with_ties = FALSE)

# Drop unused samples from tree
tree_pruned <- drop.tip(tree1, setdiff(tree1$tip.label, best_samples$CloneA))

# Make a rename map from CloneA to Group
rename_map <- best_samples %>%
  select(CloneA, Group) %>%
  deframe()

# Relabel tree tips
tree_pruned$tip.label <- rename_map[tree_pruned$tip.label]

# Optional: remove duplicate tip labels if any remain (safety)
tree_pruned <- keep.tip(tree_pruned, unique(tree_pruned$tip.label))

tree1 <- tree_pruned


# -----------------------------
# 2. Prune to common tips
# -----------------------------
common_tips <- intersect(tree1$tip.label, tree2$tip.label)

tree1_pruned <- drop.tip(tree1, setdiff(tree1$tip.label, common_tips))
tree2_pruned <- drop.tip(tree2, setdiff(tree2$tip.label, common_tips))




# # -----------------------------
# # 3. Reorder tips to match
# # --------
# library(ape)


# # prune to common tips
# common_tips <- intersect(tree1$tip.label, tree2$tip.label)
# tree1_pruned <- drop.tip(tree1, setdiff(tree1$tip.label, common_tips))
# tree2_pruned <- drop.tip(tree2, setdiff(tree2$tip.label, common_tips))


# df1 <- read.csv("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/genomic_types.csv")
# df2 <- read.csv("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/mito_types.csv")


# genomictypes <- df1[, c(2,3)]
# mitotypes <- df2[, c(2,3)]

# types <- merge(genomictypes,mitotypes, by="CloneA")

# AC_C_Clones <- subset(types, Group.x == "AC" | Group.x == "C")$CloneA
# AC_C_Clones <- subset(types, Group.y == "B" | Group.y == "E"| Group.y == "F")$CloneA


# tree1_pruned <- drop.tip(tree1_pruned, setdiff(tree1$tip.label, AC_C_Clones))
# tree2_pruned <- drop.tip(tree2_pruned, setdiff(tree2$tip.label, AC_C_Clones))


# association: one-to-one, same labels
assoc <- cbind(tree1_pruned$tip.label, tree1_pruned$tip.label)

tree1_prunedrooted_tree <- root(tree1_pruned, outgroup = "A", resolve.root = TRUE)
tree2_prunedrooted_tree <- root(tree2_pruned, outgroup = "A", resolve.root = TRUE)


cophylo<-cophylo(tree1_prunedrooted_tree,tree2_prunedrooted_tree,
    assoc=assoc)

png("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/tree_cophylo_beast_vs_neighbor.png",
    width = 6000, height = 6000, res = 300)
plot(cophylo)
dev.off()




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

Pool <- subset(metadata_with_clone, accuratelocation == "P66")

samples_to_keep <- Pool %>% pull(Well)


library(RColorBrewer)

ngroups <- length(unique(df1$Group))

# Create an extended palette (beyond 8 colors)
group_colors <- setNames(
  colorRampPalette(brewer.pal(8, "Set2"))(ngroups),
  sort(unique(df1$Group))
)

# Add the color column to df1
df1$color <- group_colors[df1$Group]

df2$color <- group_colors[df2$Group]  # or a different palette if you want










# ----------------------------
# Tip colors
# ----------------------------
cols1 <- setNames(df1$color, df1$CloneA)
cols2 <- setNames(df2$color, df2$CloneA)

tipcols1 <- cols1[match(cophylo$trees[[1]]$tip.label, names(cols1))]
tipcols2 <- cols2[match(cophylo$trees[[2]]$tip.label, names(cols2))]

tipcols1[is.na(tipcols1)] <- "black"
tipcols2[is.na(tipcols2)] <- "black"


library(phangorn) 
# Make sure your tip colors are ready
cols1 <- setNames(df1$color, df1$CloneA)
cols2 <- setNames(df2$color, df2$CloneA)

tipcols1 <- cols1[match(cophylo$trees[[1]]$tip.label, names(cols1))]
tipcols2 <- cols2[match(cophylo$trees[[2]]$tip.label, names(cols2))]

tipcols1[is.na(tipcols1)] <- "black"
tipcols2[is.na(tipcols2)] <- "black"

# Define the function
get_edge_colors <- function(tree, tipcols) {
  edge_colors <- rep("black", nrow(tree$edge))
  
  for (i in seq_len(nrow(tree$edge))) {
    node <- tree$edge[i, 2]
    # Use phangorn::Descendants for tips
    tips <- phangorn::Descendants(tree, node, "tips")[[1]]
    edge_colors[i] <- tipcols[tips[1]]  # use first descendant tip's color
  }
  return(edge_colors)
}

# Generate edge colors
edgecols1 <- get_edge_colors(cophylo$trees[[1]], tipcols1)
edgecols2 <- get_edge_colors(cophylo$trees[[2]], tipcols2)

# Combine into edge.col for cophylo
edge.col <- list(left = edgecols1, right = edgecols2)



genomictypes <- df1[, c(2,3)]
mitotypes <- df2[, c(2,3)]

types <- merge(genomictypes,mitotypes, by="CloneA")

types$types <- paste(types$Group.x, types$Group.y, sep ="_")

table(types$types)
table(types$Group.x)
table(types$Group.y)

type_counts <- as.data.frame(table(types$types))

# Rename columns
colnames(type_counts) <- c("Type", "Count")

# View the table
print(type_counts)

# Optional: sort by count descending
type_counts <- type_counts[order(-type_counts$Count), ]


# Get unique group-color mappings for left and right
left_colors  <- unique(df1[, c("Group", "color")])
right_colors <- unique(df2[, c("Group", "color")])

#left_colors <- subset(left_colors, Group == "C" | Group == "AC")
#right_colors <- subset(right_colors, Group == "B" | Group == "E"| Group == "F")

cb_palette <- data.table(
  "A" = "#F0E442",  # yellow
  "B" = "#E69F00",  # orange
  "C" = "#0072B2",  # dark blue
  "D" = "#D55E00",  # reddish-orange
  "E" = "#56B4E9",  # sky blue
  "F" = "#009E73"  # bluish-green  
)

cb_df <- data.frame(
  Group = names(cb_palette),
  color = as.character(cb_palette),
  row.names = NULL
)


# Make sure they are sorted in the order you like
left_colors  <- left_colors[order(left_colors$Group), ]
right_colors <- right_colors[order(right_colors$Group), ]



old_to_group <- setNames(right_colors$color, right_colors$Group)

# Apply it to cb_df by matching the Group
cb_df$new_color <- old_to_group[cb_df$Group]

old_to_new <- setNames(cb_df$color,cb_df$new_color)

# Replace colors in edge.col$right
edge.col$right <- old_to_new[edge.col$right]


# Check
head(edge.col$right, 20)






# Suppose cb_palette is:
# B "#E69F00"  E "#56B4E9"  F "#009E73"



metadata_rockpools <- subset(metadata_with_clone, accuratelocation == "P58"|accuratelocation == "P62"|
accuratelocation =="P63"|accuratelocation == "P66"|accuratelocation == "Gilmer")

metadata_rockpools2 <- metadata_rockpools[, c(2,9)]

# Define pond palette for your 5 sites
pond_palette <- c(
  "P58"   = "#E69F00",  # orange
  "P62"   = "#56B4E9",  # sky blue
  "P63"   = "#009E73",  # bluish-green
  "P66"   = "#D55E00",  # vermillion red
  "Gilmer"= "#CC79A7"   # purple/magenta
)

left_tips  <- cophylo$trees[[1]]$tip.label

# Make a lookup: Well -> pond
well_to_pond <- setNames(metadata_rockpools2$accuratelocation,
                         metadata_rockpools2$Well)

# Get pond for each tip in the tree
ltip_colors <- pond_palette[well_to_pond[left_tips]]


# Create a named vector: old tip label -> new label
tip_rename <- setNames(metadata_rockpools2$accuratelocation, metadata_rockpools2$Well)

png("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/cophylo_MitoBEF.png",
    res = 300, width = 4000, height = 5000)

plot(cophylo, edge.col = edge.col, fsize = 0.8, lwd = 2, link.col = "gray")

# Add legend for left tree (top left)
legend("topleft",
       legend = left_colors$Group,
       col = left_colors$color,
       lwd = 2,
       cex = 1.2,
       bty = "n",
       title = "Left tree groups")

# Add legend for right tree (top right)
legend("topright",
       legend = cb_df$Group,
       col = cb_df$color,
       lwd = 2,
       cex = 1.2,
       bty = "n",
       title = "Right tree groups")

dev.off()





# Keep original tip labels
left_tips  <- cophylo$trees[[1]]$tip.label
right_tips <- cophylo$trees[[2]]$tip.label

# Map tip labels to pond colors
well_to_pond <- setNames(metadata_rockpools2$accuratelocation,
                         metadata_rockpools2$Well)

ltip_colors <- pond_palette[well_to_pond[left_tips]]
rtip_colors <- pond_palette[well_to_pond[right_tips]]







left_tips  <- cophylo$trees[[1]]$tip.label
right_tips  <- cophylo$trees[[2]]$tip.label

# Make a lookup: Well -> pond
well_to_pond <- setNames(metadata_rockpools2$accuratelocation,
                         metadata_rockpools2$Well)

# Get pond for each tip in the tree
welltopond2 <- well_to_pond[left_tips]
welltopond3 <- well_to_pond[right_tips]



# Plot cophylo with original tip labels
png("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/cophylo_colors_BEF.png",
    res = 300, width = 4000, height = 5000)

plot(cophylo, edge.col = edge.col, fsize = 0.8, lwd = 2, link.col = "gray60", ftype="i")

# Add colored tip labels
tiplabels.cophylo(
  text = welltopond2,
  col = ltip_colors,
  frame = "none",
    adj = -1,
  which = "left"
)

tiplabels.cophylo(
  text = welltopond3,
  col = rtip_colors,
  frame = "none",
    adj = 1.5 ,
  which = "right"
)

legend("top",                     # position, adjust as needed
       legend = unique(names(rtip_colors)),   # pond names
       col = unique(rtip_colors),             # matching colors
       pch = 19,                       # point type to match tip symbols
       cex = 1.5,                      # text/point size
       bty = "n")                      # no box around legend

dev.off()


