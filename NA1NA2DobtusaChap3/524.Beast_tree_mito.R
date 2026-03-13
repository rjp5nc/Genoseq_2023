#ijob -A berglandlab -c10 -p standard --mem=40G
#module load gcc/11.4.0  openmpi/4.1.4 icu R/4.3.1;R
#R

library(dplyr)
library(ggtree)
library(treeio)
library(ape)
library(stringr)
library(ggplot2)

options(ignore.negative.edge = TRUE)

## ------------------------
## Inputs
## ------------------------
metadata_file  <- "/project/berglandlab/Robert/UKSequencing2022_2024/old_stuff/metadata_with_clone.csv"
samples_file   <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/Sample_ID_Species_merged_20251227.csv"
mitotypes_file <- "/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/mito_types.csv"

mcc_file     <- "/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/gatk_gvcf_top3_plusA/snapp.mito.highest1.tree"
out_png_fan  <- "/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/gatk_gvcf_top3_plusA/usobtusa_mito_tree_beastcirc.png"
out_pdf_fan  <- "/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/gatk_gvcf_top3_plusA/usobtusa_mito_tree_beastcirc.pdf"
outgroup_tip <- "Gilmer5_F2_clone"

## ------------------------
## Read MCC tree WITH annotations
## ------------------------
td <- read.beast(mcc_file)   # treedata keeps posterior/HPD annotations

## Root phylo inside td (keeps annotations in td)
phy <- td@phylo
if (outgroup_tip %in% phy$tip.label) {
  phy <- root(phy, outgroup = outgroup_tip, resolve.root = TRUE)
}
if (!is.null(phy$edge.length)) phy$edge.length[phy$edge.length < 0] <- 0
td@phylo <- phy

## ------------------------
## Metadata mappings
## ------------------------
metadata <- read.csv(metadata_file) %>%
  filter(!clone %in% c("Blank", "BLANK")) %>%
  mutate(
    Well = trimws(Well),
    clone = paste0(Well, "_clone"),
    accuratelocation = trimws(accuratelocation)
  ) %>%
  select(clone, accuratelocation)

samples2 <- read.csv(samples_file) %>%
  mutate(
    Sample_ID     = trimws(Sample_ID),
    Sample_ID_old = trimws(Sample_ID_old),
    clone = ifelse(str_detect(Sample_ID, "_clone$"),
                   Sample_ID,
                   paste0(Sample_ID, "_clone"))
  ) %>%
  select(clone, Sample_ID_old)

mitotypes <- read.csv(mitotypes_file) %>%
  mutate(
    sampleA = trimws(sampleA),
    Group   = trimws(Group),
    clone   = paste0(sampleA, "_clone")
  ) %>%
  select(clone, mitotype = Group)

## ------------------------
## Tip table: label + mitotype
## Label rule: Sample_ID_old > accuratelocation > clone
## ------------------------
tip_df <- tibble(clone = td@phylo$tip.label) %>%
  left_join(metadata,  by = "clone") %>%
  left_join(samples2,  by = "clone") %>%
  left_join(mitotypes, by = "clone") %>%
  mutate(
    Sample_ID_old    = na_if(Sample_ID_old, ""),
    accuratelocation = na_if(accuratelocation, ""),
    tip_label = coalesce(accuratelocation, Sample_ID_old, clone),
    mitotype  = ifelse(is.na(mitotype) | mitotype == "", "Unknown", mitotype)
  )

## ------------------------
## Plot (circular tree, colored by mitotype)
## ------------------------
p0 <- suppressWarnings(
  ggtree(td, layout = "circular") %<+% tip_df
)


node_df <- p0$data %>%
  dplyr::filter(!isTip, !is.na(posterior))

p <- p0 +
  geom_tiplab(aes(label = tip_label, color = mitotype), size = 6, offset = 0.005) +
  theme_tree2() +
  labs(color = "Mitotype") +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 16),
    legend.text  = element_text(size = 14)
  )+
  geom_text(
    data = node_df,
    aes(x = x, y = y, label = sprintf("%.2f", posterior)),
    size  = 5.0,   # large node support text
    vjust = -0.3
  )



## ------------------------
## Save
## ------------------------
ggsave(out_png_fan, plot = p, width = 14, height = 14, dpi = 300)
ggsave(out_pdf_fan, plot = p, width = 14, height = 14, dpi = 300)
