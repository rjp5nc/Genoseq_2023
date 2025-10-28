#ijob -A berglandlab -c10 -p standard --mem=40G
#module load gcc/11.4.0  openmpi/4.1.4 icu R/4.3.1
#R

library(SeqArray)
library(SNPRelate)
library(ape)
library(ggtree)
library(ggplot2)
library(dplyr)
library(tibble)
library(data.table)
library(igraph)

samples <- read.csv("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/Sample_ID_Species_merged_20250627.csv")
usobtusasamps <- subset(samples, Species == "Daphnia obtusa" & Continent == "NorthAmerica")
metadata_with_clone <- read.csv("/project/berglandlab/Robert/UKSequencing2022_2024/old_stuff/metadata_with_clone.csv", header = TRUE)
metadata_with_clone <- subset(metadata_with_clone, clone !="Blank" & clone !="BLANK")
mitotypes <- read.csv("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/mito_types.csv")


# 1. Load a single tree (like your MCC summary tree)


tree <- read.nexus("/scratch/rjp5nc/snapp5/snapp.mito.highest1.tree")
rooted_tree <- root(tree, outgroup = "Gilmer5_F2_clone", resolve.root = TRUE)

# --- 2. Prepare metadata ---
metadata_with_clone <- metadata_with_clone %>%
  mutate(clone = paste0(Well, "_clone")) %>%
  filter(location %in% c("Rockpool", "US"))

# --- 3. Create clone → Well mapping ---
label_map <- metadata_with_clone %>%
  select(clone, Well) %>%
  distinct() %>%
  mutate(
    clone = trimws(as.character(clone)),
    Well = trimws(as.character(Well))
  ) %>%
  deframe()

# --- 4. Relabel tree tips with Well IDs ---
rooted_tree$tip.label <- trimws(as.character(rooted_tree$tip.label))
rooted_tree$tip.label <- label_map[rooted_tree$tip.label]

# Check for unmatched tips
missing <- sum(is.na(rooted_tree$tip.label))
if (missing > 0) {
  warning(paste(missing, "tree tips did not match any 'clone' in metadata_with_clone"))
}

# --- 5. Prepare group color info ---
mitotypes <- mitotypes %>%
  mutate(clone_paste = paste0(CloneA, "_clone"))

tip_color_dt <- data.table(
  Group = mitotypes$Group,
  clone = mitotypes$clone_paste
)

# Map clone names to Well names
tip_color_dt <- tip_color_dt %>%
  mutate(Well = label_map[clone]) %>%
  filter(!is.na(Well))

cb_palette <- c(
  "#E69F00", # orange
  "#56B4E9", # sky blue
  "#009E73", # green
  "#6e6923", # yellow
  "#034063", # blue
  "#D55E00", # red-orange
  "#CC79A7", # pink
  "#000000"  # black
)



# Get unique Groups
groups <- unique(tip_color_dt$Group)

# Assign colors to groups (recycle palette if more groups than colors)
group_colors <- setNames(cb_palette[1:length(groups)], groups)


# Create a named color vector for the tree tips
tip_colors <- setNames(group_colors[tip_color_dt$Group], tip_color_dt$Well)
tip_colors <- tip_colors[rooted_tree$tip.label]  # reorder to match tip order



label_map <- metadata_with_clone %>%
  select(clone, accuratelocation) %>%
  distinct() %>%
  mutate(
    clone = trimws(as.character(clone)),
    accuratelocation = trimws(as.character(accuratelocation))
  ) %>%
  deframe()  # names = clone, values = accuratelocation

# --- 4. Relabel tree tips using accurate_location ---
rooted_tree$tip.label <- trimws(as.character(rooted_tree$tip.label))

label_map2 <- setNames(label_map, gsub("_clone$", "", names(label_map)))

# Then replace tip labels
rooted_tree$tip.label <- label_map2[rooted_tree$tip.label]

# --- 6. Save a PDF of the tree (simple version) ---
pdf("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/beast_mito_tree_MCC.pdf",
    width = 20, height = 40)
plot(rooted_tree, main = "MCC Tree")
dev.off()

# --- 7. Save a colored circular PNG plot ---
png("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/usobtusa_mito_tree_beastcirc.png",
    res = 300, width = 2000, height = 2000)

plot.phylo(rooted_tree,
           type = "fan",
           cex = 0.8,
           label.offset = 0.01,
           no.margin = TRUE,
           tip.color = tip_colors,
           main = "Rooted Mito Tree (root = Gilmer5_F2_clone)")

legend("topleft",
       legend = names(group_colors),
       col = group_colors,
       pch = 19,
       pt.cex = 1.5,
       cex = 1,
       bty = "n",
       title = "Group")

dev.off()






node_ages_subs <- branching.times(rooted_tree)

# convert to absolute time
mu <- 1.75e-7          # mutation rate per site per generation
gens_per_year <- 6     # adjust to your species if known

node_ages_gen  <- node_ages_subs / mu
node_ages_year <- node_ages_gen / gens_per_year

# combine into a table
node_age_df <- data.frame(
  node = names(node_ages_subs),
  subs_per_site = node_ages_subs,
  generations   = node_ages_gen,
  years         = node_ages_year
)

head(node_age_df)



is.ultrametric(rooted_tree)






mu <- 1.75e-7                 # mutation rate per site per generation
gens_per_year <- 6           

# === Read tree ===

tree2 <- tree


tree$edge.length
# === Convert branch lengths from substitutions/site to generations ===
tree2$edge.length <- tree2$edge.length / mu

# Get node ages in generations
node_ages_gen <- branching.times(tree2)

# Convert to years
node_ages_years <- node_ages_gen / gens_per_year

# Combine into a data frame
node_age_df <- data.frame(
  node = names(node_ages_gen),
  generations = as.numeric(node_ages_gen),
  years = as.numeric(node_ages_years)
) %>%
  arrange(desc(generations))

# View top nodes
head(node_age_df)

node_age_df$yearsnormalized <- node_age_df$years / 6542.068




png("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/usobtusa_mito_tree_beastcirc_years.png",
    res = 300, width = 2000, height = 2000)

plot.phylo(rooted_tree,
           type = "fan",
           cex = 0.8,
           label.offset = 0.01,
           no.margin = TRUE,
           tip.color = tip_colors,
           main = "Rooted Mito Tree (root = Gilmer5_F2_clone)")

# --- Add internal node labels: years ---
# ape internal nodes are numbered Ntip+1:Nnode
Ntip <- length(rooted_tree$tip.label)
for(i in 1:nrow(node_age_df)){
  node_num <- as.numeric(node_age_df$node[i])
  # Only label internal nodes
  if(node_num > Ntip){
    nodelabels(text = round(node_age_df$yearsnormalized[i], 0), 
               node = node_num, 
               frame = "none", 
               adj = c(0.5, -0.5), 
               cex = 0.6, 
               col = "darkred")
  }
}

legend("topleft",
       legend = names(group_colors),
       col = group_colors,
       pch = 19,
       pt.cex = 1.5,
       cex = 1,
       bty = "n",
       title = "Group")

dev.off()
