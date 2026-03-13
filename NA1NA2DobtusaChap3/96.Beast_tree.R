#ijob -A berglandlab -c10 -p standard --mem=40G
#module load gcc/11.4.0  openmpi/4.1.4 icu R/4.3.1
#R

library(ape)

# 1. Load a single tree (like your MCC summary tree)

tree <- read.nexus("/scratch/rjp5nc/snapp5/snapp2.tree")
rooted_tree <- root(tree, outgroup = "A", resolve.root = TRUE)


pdf("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/beast_tree_MCC.pdf",
    width = 20, height = 40)

plot(rooted_tree)  # quick visualization
title("MCC Tree")

dev.off()





# 2. Load a set of trees (posterior sample)
trees <- read.nexus("/scratch/rjp5nc/snapp5/snapp.mono2.trees")
length(trees)   # number of posterior trees


cons <- consensus(trees, p = 0.5)  # 50% majority rule


pdf("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/beast_cons_tree.pdf",
    width = 15, height = 30)
plot(cons)
nodelabels(
  text = round(cons$node.label, 2),  # node support
  cex = 1.5,                          # label size
  frame = "rect",                      # rectangle box around label
  adj = c(0.5, 0),                     # 0 = above node, 0.5 = horizontally centered
  font = 2                             # bold text (optional)
)
dev.off()



node_ages <- max(node.depth.edgelength(cons)) - node.depth.edgelength(cons)


n_tips <- Ntip(cons)
internal_nodes <- (n_tips + 1):(n_tips + cons$Nnode)

# Create a table of internal node ages
node_table <- data.frame(
  node = internal_nodes,
  age = node_ages[internal_nodes]
)
head(node_table)





pdf("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/beast_diverge_tree.pdf",
    width = 20, height = 40)

plot(cons, cex = 1)  # tip labels bigger if needed

# Add node ages in boxes over nodes
nodelabels(
  text = round(node_table$age, 3),
  node = node_table$node,
  frame = "rect",
  bg = "lightblue",
  cex = 1.2,
  adj = c(0.5, 0)
)

title("SNAPP Tree with Node Divergence Times")
dev.off()








tree <- trees[[1]]
n_tips <- Ntip(tree)
n_nodes <- tree$Nnode

# Compute node heights (distance from root to each node)
# Use nodeHeights from phangorn package for rooted trees
library(phangorn)
library(phytools)

heights <- nodeHeights(tree)  # matrix: start and end height of each edge
# Node age = max tree height - node height at top of edge
tree_height <- max(heights)
node_ages <- numeric(n_nodes)



# Example: tree already loaded
n_tips <- Ntip(tree)
n_nodes <- tree$Nnode

# Function to compute node heights (distance from root)
nodeHeights2 <- function(tree) {
  eh <- matrix(NA, nrow = Nedge(tree <- tree)$edge[,1] %>% length(), ncol = 2)
  h <- numeric(n_tips + n_nodes) # node heights
  for(i in 1:nrow(tree$edge)){
    parent <- tree$edge[i,1]
    child  <- tree$edge[i,2]
    start <- h[parent]
    end   <- start + tree$edge.length[i]
    eh[i,] <- c(start, end)
    h[child] <- end
  }
  list(edgeHeights = eh, nodeHeights = h)
}

res <- nodeHeights2(tree)
h_nodes <- res$nodeHeights
tree_height <- max(h_nodes)

# Internal nodes
internal_nodes <- (n_tips+1):(n_tips+n_nodes)
node_ages <- tree_height - h_nodes[internal_nodes]

node_table <- data.frame(
  node = internal_nodes,
  age  = node_ages
)
head(node_table)







pdf("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/beast_diverge_tree.pdf",
    width = 20, height = 40)

plot(tree, cex = 1)  # increase cex for bigger tip labels if needed

# Add node ages in boxes above internal nodes
nodelabels(
  text = round(node_table$age, 5),
  node = node_table$node,
  frame = "rect",      # rectangle box
  bg = "lightblue",    # box color
  cex = 1.2,           # font size
  adj = c(0.5, 0)      # center horizontally, put above node
)

title("SNAPP Tree with Node Divergence Times")

dev.off()



mu <- 3.8e-9        # substitutions per site per generation
gen_time <- 0.20      # years per generation

# Convert node ages to years
node_table$age_years <- node_table$age / mu * gen_time

# Check
head(node_table)

rooted_tree <- root(tree, outgroup = "T", resolve.root = TRUE)


pdf("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/beast_diverge_tree.pdf",
    width = 20, height = 40)

plot(rooted_tree, cex = 1)  # increase cex for bigger tip labels if needed

# Add node ages in boxes above internal nodes
nodelabels(
  text = round(node_table$age_years, 1),
  node = node_table$node,
  frame = "rect",      # rectangle box
  bg = "lightblue",    # box color
  cex = 1.2,           # font size
  adj = c(0.5, 0)      # center horizontally, put above node
)

title("SNAPP Tree with Node Divergence Times")

dev.off()




library(tidyverse)

# Keep only valid phylo objects
trees <- trees[sapply(trees, function(x) inherits(x, "phylo"))]

n_trees <- length(trees)
n_tips <- Ntip(trees[[1]])
n_nodes <- trees[[1]]$Nnode
internal_nodes <- (n_tips+1):(n_tips+n_nodes)

all_node_ages <- matrix(NA, nrow = n_nodes, ncol = n_trees)

nodeHeights2 <- function(tree){
  if(!inherits(tree, "phylo")) stop("tree must be phylo")
  n <- Ntip(tree) + tree$Nnode
  h <- numeric(n)               # node heights
  for(i in 1:nrow(tree$edge)){
    parent <- tree$edge[i,1]
    child  <- tree$edge[i,2]
    start <- h[parent]
    end   <- start + tree$edge.length[i]
    h[child] <- end
  }
  h
}

for(t in 1:n_trees){
  h_nodes <- nodeHeights2(trees[[t]])
  tree_height <- max(h_nodes)
  all_node_ages[,t] <- tree_height - h_nodes[internal_nodes]
}

mean_ages <- rowMeans(all_node_ages, na.rm = TRUE)

# 95% HPD (using quantiles)
HPD_lower <- apply(all_node_ages, 1, function(x) quantile(x, 0.025, na.rm = TRUE))
HPD_upper <- apply(all_node_ages, 1, function(x) quantile(x, 0.975, na.rm = TRUE))

node_table <- data.frame(
  node = internal_nodes,
  mean_age = mean_ages,
  HPD_lower = HPD_lower,
  HPD_upper = HPD_upper
)

head(node_table)


mu <- 3.8e-9
gen_time <- 0.25

node_table$mean_age_years <- node_table$mean_age / mu * gen_time
node_table$HPD_lower_years <- node_table$HPD_lower / mu * gen_time
node_table$HPD_upper_years <- node_table$HPD_upper / mu * gen_time