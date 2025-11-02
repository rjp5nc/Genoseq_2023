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
tree2_file <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/usdobtusa_tree_mito.nwk"

tree1 <- read.tree(tree1_file)
tree2 <- read.tree(tree2_file)

# -----------------------------
# 2. Prune to common tips
# -----------------------------
common_tips <- intersect(tree1$tip.label, tree2$tip.label)

tree1_pruned <- drop.tip(tree1, setdiff(tree1$tip.label, common_tips))
tree2_pruned <- drop.tip(tree2, setdiff(tree2$tip.label, common_tips))

# -----------------------------
# 3. Reorder tips to match
# --------
library(ape)


# prune to common tips
common_tips <- intersect(tree1$tip.label, tree2$tip.label)
tree1_pruned <- drop.tip(tree1, setdiff(tree1$tip.label, common_tips))
tree2_pruned <- drop.tip(tree2, setdiff(tree2$tip.label, common_tips))

# association: one-to-one, same labels
assoc <- cbind(tree1_pruned$tip.label, tree1_pruned$tip.label)

# save bigger PNG
png("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/tree_cophylo_bold_only2.png",
    width = 6000, height = 12000, res = 300)

cophyloplot(tree1_pruned, tree2_pruned,
            assoc = assoc,
            use.edge.length = FALSE,
            space = 50,                # less gap → bigger trees
            gap = 200,                   # closer to the axis
            length.line = 100,           # connecting line length
            col = "gray40",            # connector line color
            lwd = 1,                   # thicker tree branches
            cex = 1.2,                 # tip label size
            direction = c("right","left"))  # mirrored orientation

dev.off()


cophylo<-cophylo(tree1_pruned,tree2_pruned,
    assoc=assoc)

png("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/tree_cophylo_bold_only2.png",
    width = 6000, height = 12000, res = 300)
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


df1 <- read.csv("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/genomic_types.csv")
df2 <- read.csv("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/mito_types.csv")

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
right_tips <- cophylo$trees[[2]]$tip.label

mito_pallete <- c(
  "A"   = "#66C2A5", 
  "B"   = "#E59077",
  "C"   = "#9C9CC9",  
  "D"   = "#AF9AAA", 
  "E"   = "#CA9591",   
  "F"   = "#949EC3"  
)


# Make a lookup: Well -> pond
well_to_pond <- setNames(metadata_rockpools2$accuratelocation,
                         metadata_rockpools2$Well)

# Get pond for each tip in the tree
ltip_colors <- pond_palette[well_to_pond[left_tips]]
rtip_colors <- pond_palette[well_to_pond[right_tips]]


png("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/cophylo_colors_only2.png", res = 300, width = 4000, height = 19000)

plot(cophylo, edge.col = edge.col, fsize = 0.8, lwd = 2, link.col = ltip_colors)

dev.off()


# Save to CSV

library(tidyverse)


left_tips  <- cophylo$trees[[1]]$tip.label
right_tips  <- cophylo$trees[[2]]$tip.label

# Make a lookup: Well -> pond
well_to_pond <- setNames(metadata_rockpools2$accuratelocation,
                         metadata_rockpools2$Well)

# Get pond for each tip in the tree
welltopond2 <- well_to_pond[left_tips]
welltopond3 <- well_to_pond[right_tips]



png("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/cophylo_colors.png",
    res = 1000, width = 7000, height = 8000)

# Plot cophylo with colored edges
plot(cophylo, edge.col = edge.col, fsize = 0.8, lwd = 2, link.col = "gray40")

# Add legend for left tree tip colors
legend("top",                     # position, adjust as needed
       legend = names(pond_palette),   # pond names
       col = pond_palette,             # matching colors
       pch = 19,                       # point type to match tip symbols
       cex = 1.5,                      # text/point size
       bty = "n")                      # no box around legend
                   # no box around legend

legend("topright",                     # position, adjust as needed
       legend = names(mito_pallete),   # pond names
       col = mito_pallete,             # matching colors
       pch = 19,                       # point type to match tip symbols
       cex = 1.5,                      # text/point size
       bty = "n")

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


dev.off()



