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
mitotypes <- read.csv("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/mito_types_v2.csv")




# 1. Load a single tree (like your MCC summary tree)


tree <- read.nexus("/scratch/rjp5nc/snapp5/snapp.mito.highest1.tree")
rooted_tree <- root(tree, outgroup = "Gilmer5_F2_clone", resolve.root = TRUE)






library(treeio)
library(ggtree)

mcc_file <- "/scratch/rjp5nc/snapp5/snapp.mito.highest1.tree"

tr <- read.beast(mcc_file)   # keeps posterior etc.

# circular plot with posterior on internal nodes
p <- ggtree(tr, layout = "fan") %<+% tr@data +
  geom_tiplab(size = 2) +
  geom_text2(aes(subset = !isTip, label = sprintf("%.2f", posterior)),
             size = 2, vjust = -0.3)

p





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

mito_wells <- tip_color_dt$Well
mito_wells <- unique(mito_wells)

tips_to_drop <- setdiff(rooted_tree$tip.label, mito_wells)

tree_mito_only <- drop.tip(rooted_tree, tips_to_drop)



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
tip_colors <- tip_colors[tree_mito_only$tip.label]  # reorder to match tip order



label_map <- metadata_with_clone %>%
  select(clone, accuratelocation) %>%
  distinct() %>%
  mutate(
    clone = trimws(as.character(clone)),
    accuratelocation = trimws(as.character(accuratelocation))
  ) %>%
  deframe()  # names = clone, values = accuratelocation

# --- 4. Relabel tree tips using accurate_location ---
tree_mito_only$tip.label <- trimws(as.character(tree_mito_only$tip.label))


label_map2 <- setNames(label_map, gsub("_clone$", "", names(label_map)))


# Then replace tip labels
tree_mito_only$tip.label <- label_map2[tree_mito_only$tip.label]


# --- 6. Save a PDF of the tree (simple version) ---
pdf("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/beast_mito_tree_MCC.pdf",
    width = 20, height = 40)
plot(tree_mito_only, main = "MCC Tree")
dev.off()

# --- 7. Save a colored circular PNG plot ---
png("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/usobtusa_mito_tree_beastcirc.png",
    res = 300, width = 2000, height = 2000)

plot.phylo(tree_mito_only,
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

















library(dplyr)
library(data.table)
library(ape)
library(treeio)
library(ggtree)

## ------------------------
## 0) Inputs
## ------------------------
samples_file   <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/Sample_ID_Species_merged_20250627.csv"
metadata_file  <- "/project/berglandlab/Robert/UKSequencing2022_2024/old_stuff/metadata_with_clone.csv"
mitotypes_file <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/mito_types_v2.csv"
mcc_file       <- "/scratch/rjp5nc/snapp5/snapp.mito.highest1.tree"

out_pdf_rect   <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/beast_mito_tree_MCC.pdf"
out_png_fan    <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/usobtusa_mito_tree_beastcirc.png"

outgroup_tip   <- "Gilmer5_F2_clone"   # outgroup used for rooting (in original clone naming)

## ------------------------
## 1) Load data
## ------------------------
samples <- read.csv(samples_file)
usobtusasamps <- subset(samples, Species == "Daphnia obtusa" & Continent == "NorthAmerica")

metadata_with_clone <- read.csv(metadata_file, header = TRUE)
metadata_with_clone <- subset(metadata_with_clone, clone != "Blank" & clone != "BLANK")

mitotypes <- read.csv(mitotypes_file)

## ------------------------
## 2) Read MCC tree (BEAST-aware) and root
## ------------------------
tr <- read.beast(mcc_file)  # treedata, keeps posterior etc.

# Root the phylo inside treedata
tr@phylo <- root(tr@phylo, outgroup = outgroup_tip, resolve.root = TRUE)

## ------------------------
## 3) Prepare metadata and maps
## ------------------------
# Add clone name used in the tree pipeline
metadata_with_clone <- metadata_with_clone %>%
  mutate(clone = paste0(Well, "_clone")) %>%
  filter(location %in% c("Rockpool", "US")) %>%
  mutate(
    clone = trimws(as.character(clone)),
    Well  = trimws(as.character(Well)),
    accuratelocation = trimws(as.character(accuratelocation))
  )

# clone -> Well map
clone_to_well <- metadata_with_clone %>%
  select(clone, Well) %>%
  distinct() %>%
  deframe()

# clone -> accuratelocation map (and build a version without "_clone" suffix)
clone_to_loc <- metadata_with_clone %>%
  select(clone, accuratelocation) %>%
  distinct() %>%
  deframe()

clone_to_loc2 <- setNames(clone_to_loc, gsub("_clone$", "", names(clone_to_loc)))

## ------------------------
## 4) Prepare mitotype → tip colors
## ------------------------
mitotypes <- mitotypes %>%
  mutate(clone_paste = paste0(CloneA, "_clone"))

tip_color_dt <- data.frame(
  Group = mitotypes$Group,
  clone = mitotypes$clone_paste,
  stringsAsFactors = FALSE
) %>%
  mutate(
    clone = trimws(as.character(clone)),
    Well  = clone_to_well[clone]
  ) %>%
  filter(!is.na(Well)) %>%
  distinct(Group, Well)

# Color-blind friendly palette (recycled if needed)
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

groups <- unique(tip_color_dt$Group)
group_colors <- setNames(cb_palette[seq_len(min(length(cb_palette), length(groups)))], groups)

# If more groups than palette length, recycle safely
if (length(groups) > length(cb_palette)) {
  group_colors <- setNames(rep(cb_palette, length.out = length(groups)), groups)
}

# named vector: names are Well, values are hex colors
tip_colors_by_well <- setNames(group_colors[tip_color_dt$Group], tip_color_dt$Well)

# wells that have mitotypes (used for pruning)
mito_wells <- unique(tip_color_dt$Well)

## ------------------------
## 5) Relabel tips (clone -> Well) and prune to mitotyped wells
## ------------------------
# Work on a copy of the phylo for relabeling
phy <- tr@phylo
phy$tip.label <- trimws(as.character(phy$tip.label))

# Replace clone labels with Well IDs
phy$tip.label <- clone_to_well[phy$tip.label]

# Drop any unmatched tips (NA after mapping)
na_tips <- which(is.na(phy$tip.label))
if (length(na_tips) > 0) {
  warning(paste(length(na_tips), "tips did not match any clone in metadata_with_clone; dropping them"))
  phy <- drop.tip(phy, na_tips)
}

# Now prune to only wells found in mitotypes
tips_to_drop <- setdiff(phy$tip.label, mito_wells)
tree_mito_only <- drop.tip(phy, tips_to_drop)

# Colors aligned to the pruned tree
tip_colors_mito <- tip_colors_by_well[tree_mito_only$tip.label]
if (any(is.na(tip_colors_mito))) {
  bad <- tree_mito_only$tip.label[is.na(tip_colors_mito)]
  stop("Some remaining tips have no color (unexpected). Examples: ", paste(head(bad), collapse = ", "))
}

## ------------------------
## 6) Make a pruned treedata that matches tree_mito_only
##    (keeps posterior in tr@data)
## ------------------------
keep_tips <- tree_mito_only$tip.label

# Rebuild tr_mito by pruning using the original treedata labels (clone form),
# then relabel to wells again to match tree_mito_only.
# First, figure out which original clone tips correspond to keep_tips wells:
clone_keep <- names(clone_to_well)[clone_to_well %in% keep_tips]
clone_keep <- trimws(as.character(clone_keep))

drop_clone <- setdiff(tr@phylo$tip.label, clone_keep)
tr_mito <- treeio::drop.tip(tr, tip = drop_clone)

# Relabel tr_mito tips to wells (same as above)
tr_mito@phylo$tip.label <- clone_to_well[trimws(as.character(tr_mito@phylo$tip.label))]

# Drop NAs if any
na_tips2 <- which(is.na(tr_mito@phylo$tip.label))
if (length(na_tips2) > 0) tr_mito@phylo <- drop.tip(tr_mito@phylo, na_tips2)

# Finally ensure the treedata phylo tip order equals keep_tips set
# (ggtree uses the phylo inside tr_mito)
if (!all(tr_mito@phylo$tip.label %in% keep_tips)) {
  stop("Pruned treedata still contains tips not in keep_tips; check mapping keys.")
}

## ------------------------
## 7) Save rectangular PDF (simple)
## ------------------------
pdf(out_pdf_rect, width = 20, height = 40)
plot(tr_mito@phylo, main = "MCC Tree (mitotyped only)")
dev.off()











options(ignore.negative.edge = TRUE)

## tip color + display label table
tip_df <- data.frame(
  label   = tr_mito@phylo$tip.label,
  tip_col = unname(tip_colors_by_well[tr_mito@phylo$tip.label]),
  stringsAsFactors = FALSE
)
tip_df$tip_lab <- unname(label_map2[tip_df$label])
tip_df$tip_lab[is.na(tip_df$tip_lab)] <- tip_df$label[is.na(tip_df$tip_lab)]

## base ggtree
p0 <- ggtree(tr_mito, layout = "fan")

## attach tip metadata (no branch coloring)
dd <- p0$data %>%
  left_join(tip_df, by = "label")

p0$data <- dd

## identify branch-length column
bl_col <- dplyr::case_when(
  "branch.length" %in% names(dd) ~ "branch.length",
  "length"        %in% names(dd) ~ "length",
  TRUE ~ NA_character_
)
if (is.na(bl_col)) stop("No branch length column found.")

## compute branch midpoints for posterior labels
dd <- dd %>%
  mutate(
    bl    = .data[[bl_col]],
    x_mid = x - 0.5 * bl
  )

## plot
p <- p0 +
  # draw tree with uniform branch color
  geom_tree(color = "grey40", linewidth = 0.6) +

  # tip labels colored by group
  geom_tiplab(
    aes(label = tip_lab, color = tip_col),
    size = 2.8,
    offset = 0.01,
    show.legend = FALSE
  ) +

  # posterior values on branches (not nodes)
  geom_text(
    data = dd %>% filter(!isTip),
    aes(x = x_mid, y = y, label = sprintf("%.2f", posterior)),
    inherit.aes = FALSE,
    size = 3.0,
    fontface = "bold"
  ) +

  scale_color_identity() +
  ggtitle("Rooted Mito Tree (mitotyped only; posterior on branches)")

ggsave(out_png_fan, plot = p, width = 7, height = 7, dpi = 300)
