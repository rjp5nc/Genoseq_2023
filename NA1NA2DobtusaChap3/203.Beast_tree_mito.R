#ijob -A berglandlab -c10 -p standard --mem=40G
#module load gcc/11.4.0  openmpi/4.1.4 icu R/4.3.1
#R
library(ape)
library(treeio)
library(ggtree)
library(ggplot2)
library(dplyr)
library(cowplot)

options(ignore.negative.edge = TRUE)

## ------------------------
## Inputs
## ------------------------
mcc_file       <- "/scratch/rjp5nc/snapp5/snapp.mito.highest1.tree"
metadata_file  <- "/project/berglandlab/Robert/UKSequencing2022_2024/old_stuff/metadata_with_clone.csv"
mitotypes_file <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/mito_types_v2.csv"
out_png_fan    <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/usobtusa_mito_tree_beastcirc.png"
outgroup_tip   <- "Gilmer5_F2_clone"

## ------------------------
## Read MCC tree (BEAST-aware)
## ------------------------
tr <- read.beast(mcc_file)
tr@phylo <- root(tr@phylo, outgroup = outgroup_tip, resolve.root = TRUE)

## ------------------------
## Metadata mappings
## ------------------------
metadata <- read.csv(metadata_file) %>%
  filter(clone != "Blank", clone != "BLANK") %>%
  mutate(
    clone = paste0(Well, "_clone"),
    clone = trimws(clone),
    Well  = trimws(Well),
    accuratelocation = trimws(accuratelocation)
  )

clone_to_well <- metadata %>%
  select(clone, Well) %>% distinct() %>% deframe()

well_to_loc <- metadata %>%
  select(Well, accuratelocation) %>% distinct() %>% deframe()

## ------------------------
## Mitotypes → colors
## ------------------------
mitotypes <- read.csv(mitotypes_file) %>%
  mutate(clone = paste0(CloneA, "_clone"))

tip_color_dt <- mitotypes %>%
  mutate(Well = clone_to_well[clone]) %>%
  filter(!is.na(Well)) %>%
  distinct(Group, Well)

cb_palette <- c(
  "#E69F00", "#56B4E9", "#009E73", "#6e6923",
  "#034063", "#D55E00", "#CC79A7", "#000000"
)

groups <- unique(tip_color_dt$Group)
group_colors <- setNames(rep(cb_palette, length.out = length(groups)), groups)

tip_colors_by_well <- setNames(group_colors[tip_color_dt$Group],
                               tip_color_dt$Well)

mito_wells <- unique(tip_color_dt$Well)

## ------------------------
## Prune tree to mitotyped individuals
## ------------------------
phy <- tr@phylo
phy$tip.label <- clone_to_well[trimws(phy$tip.label)]
phy <- drop.tip(phy, which(is.na(phy$tip.label)))
phy <- drop.tip(phy, setdiff(phy$tip.label, mito_wells))

keep_clones <- names(clone_to_well)[clone_to_well %in% phy$tip.label]
tr_mito <- treeio::drop.tip(tr, setdiff(tr@phylo$tip.label, keep_clones))
tr_mito@phylo$tip.label <- clone_to_well[tr_mito@phylo$tip.label]

## ------------------------
## Tip labels + colors
## ------------------------
tip_df <- data.frame(
  label   = tr_mito@phylo$tip.label,
  tip_col = tip_colors_by_well[tr_mito@phylo$tip.label],
  tip_lab = well_to_loc[tr_mito@phylo$tip.label],
  stringsAsFactors = FALSE
)
tip_df$tip_lab[is.na(tip_df$tip_lab)] <- tip_df$label[is.na(tip_df$tip_lab)]

## ------------------------
## Build tree plot
## ------------------------
p0 <- ggtree(tr_mito, layout = "fan")
dd <- p0$data %>% left_join(tip_df, by = "label")

bl_col <- if ("branch.length" %in% names(dd)) "branch.length" else "length"

dd <- dd %>%
  mutate(
    bl    = .data[[bl_col]],
    x_mid = x - 0.5 * bl
  )

p0$data <- dd

p <- p0 +
  geom_tree(color = "grey40", linewidth = 0.6) +
  geom_tiplab(
    aes(label = tip_lab, color = tip_col),
    size = 2.8,
    offset = 0.01,
    show.legend = FALSE
  ) +
  geom_text(
    data = dd %>% filter(!isTip),
    aes(x = x_mid, y = y, label = sprintf("%.2f", posterior)),
    inherit.aes = FALSE,
    size = 3.0,
    fontface = "bold"
  ) +
  scale_color_identity() +
  ggtitle("Mito Tree (rooted to Mitotype C)")

## ------------------------
## Legend (inside plot)
## ------------------------
legend_grob <- get_legend(
  ggplot(data.frame(Group = names(group_colors)),
         aes(x = Group, y = 1, color = Group)) +
    geom_point(size = 3) +
    scale_color_manual(values = group_colors, name = "Mitotype") +
    theme_void() +
    theme(legend.position = "right")
)

p_final <- ggdraw(p) +
  draw_grob(
    legend_grob,
    x = 0.005, y = 0.67,   # top-left inside plot
    width = 0.30, height = 0.36
  )

ggsave(out_png_fan, plot = p_final, width = 7, height = 7, dpi = 300)
