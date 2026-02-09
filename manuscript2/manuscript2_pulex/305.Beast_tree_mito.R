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
metadata_file  <- "/project/berglandlab/Robert/UKSequencing2022_2024/2022_2024seqmetadata20250811.csv"
samples_file   <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/Sample_ID_Species_merged_20251227.csv"
mitotypes_file <- "/scratch/rjp5nc/UK2022_2024/mitogvcf/gvcf/eudpulex_mito_reverse_out_all/locset_n38__beavercreek-birdhut-cafe__etc/eudpulex_mito_types.csv"

mcc_file     <- "/scratch/rjp5nc/UK2022_2024/redone_mito/eupulex/snapp_from_allsites/snapp.mito.highest1.tree"
out_png_fan  <- "/scratch/rjp5nc/UK2022_2024/redone_mito/eupulex/snapp_from_allsites/eupulex_mito_tree_beastcirc.png"
out_pdf_fan  <- "/scratch/rjp5nc/UK2022_2024/redone_mito/eupulex/snapp_from_allsites/eupulex_mito_tree_beastcirc.pdf"
outgroup_tip <- "RobertUK_C8_clone"




# RobertUK_C6_clone           Dorset            DRail
# 159   RobertUK_F6_clone           Dorset      DRail
# 160   RobertUK_A4_clone           Dorset      DRail
# 161  RobertUK_C12_clone           Dorset      DRail
# 162   RobertUK_A1_clone           Dorset      DRamps
# 163   RobertUK_G9_clone           Dorset      DTiny
# 164  RobertUK_F10_clone           Dorset      DRamps
# 165  RobertUK_F11_clone           Dorset      DRusty
# 166   RobertUK_F7_clone           Dorset      DSandy
# 167   RobertUK_F9_clone           Dorset      DRusty
# 168  RobertUK_G10_clone           Dorset      DTiny
# 169   RobertUK_E4_clone           Dorset      DTiny
# 170   RobertUK_C3_clone           Dorset      DTiny
# 171   RobertUK_D5_clone           Dorset      DTiny
# 172   RobertUK_G8_clone           Dorset      DTiny
# 173   RobertUK_E8_clone           Dorset      DTiny
# 127   RobertUK_G6_clone           Dorset      Doak
# 114   RobertUK_C9_clone           Dorset      DNorden
# 14    Gilmer5_H10_clone           Dorset      DIris
# 15    Gilmer5_H11_clone           Dorset      DBunk
# 16    Gilmer5_H12_clone           Dorset      DLily
# 17    RobertUK_D3_clone           Dorset      D8
# 18    RobertUK_C7_clone           Dorset      D8
# 19    RobertUK_H1_clone           Dorset      D8     
# 20    RobertUK_E9_clone           Dorset      D8
# 21   RobertUK_A12_clone           Dorset      DBunk
# 22    RobertUK_D4_clone           Dorset      DBunk
# 23   RobertUK_H11_clone           Dorset      DBunk
# 24   RobertUK_E11_clone           Dorset      DBunk
# 25    RobertUK_F2_clone           Dorset      DBunk
# 26    RobertUK_H2_clone           Dorset      DCat
# 27    RobertUK_B4_clone           Dorset      DCat
# 28    RobertUK_G4_clone           Dorset      DCat
# 29    RobertUK_F1_clone           Dorset      DCat
# 30    RobertUK_A5_clone           Dorset      DIris
# 31   RobertUK_H10_clone           Dorset      DIris
# 32   RobertUK_H12_clone           Dorset      DLily
# 33   RobertUK_C11_clone           Dorset      DLily
# 34   RobertUK_E10_clone           Dorset      DLily
# 35    RobertUK_F8_clone           Dorset      DMountie
# 36   RobertUK_G11_clone           Dorset      DMountie
# 37    RobertUK_F3_clone           Dorset      DMud
# 38    RobertUK_F4_clone           Dorset      DNorden


repl <- tribble(
  ~Well,                ~accuratelocation_new,

  # DRail
  "RobertUK_C6_clone",   "DRail",
  "RobertUK_F6_clone",   "DRail",
  "RobertUK_A4_clone",   "DRail",
  "RobertUK_C12_clone",  "DRail",

  # DRamps
  "RobertUK_A1_clone",   "DRamps",
  "RobertUK_F10_clone",  "DRamps",

  # DTiny
  "RobertUK_G9_clone",   "DTiny",
  "RobertUK_G10_clone",  "DTiny",
  "RobertUK_E4_clone",   "DTiny",
  "RobertUK_C3_clone",   "DTiny",
  "RobertUK_D5_clone",   "DTiny",
  "RobertUK_G8_clone",   "DTiny",
  "RobertUK_E8_clone",   "DTiny",

  # DRusty
  "RobertUK_F11_clone",  "DRusty",
  "RobertUK_F9_clone",   "DRusty",

  # DSandy
  "RobertUK_F7_clone",   "DSandy",

  # Doak
  "RobertUK_G6_clone",   "Doak",

  # DNorden
  "RobertUK_C9_clone",   "DNorden",
  "RobertUK_F4_clone",   "DNorden",

  # DIris
  "Gilmer5_H10_clone",   "DIris",
  "RobertUK_A5_clone",   "DIris",
  "RobertUK_H10_clone",  "DIris",

  # DBunk
  "Gilmer5_H11_clone",   "DBunk",
  "RobertUK_A12_clone",  "DBunk",
  "RobertUK_D4_clone",   "DBunk",
  "RobertUK_H11_clone",  "DBunk",
  "RobertUK_E11_clone",  "DBunk",
  "RobertUK_F2_clone",   "DBunk",

  # DLily
  "Gilmer5_H12_clone",   "DLily",
  "RobertUK_H12_clone",  "DLily",
  "RobertUK_C11_clone",  "DLily",
  "RobertUK_E10_clone",  "DLily",

  # D8
  "RobertUK_D3_clone",   "D8",
  "RobertUK_C7_clone",   "D8",
  "RobertUK_H1_clone",   "D8",
  "RobertUK_E9_clone",   "D8",

  # DCat
  "RobertUK_H2_clone",   "DCat",
  "RobertUK_B4_clone",   "DCat",
  "RobertUK_G4_clone",   "DCat",
  "RobertUK_F1_clone",   "DCat",

  # DMountie
  "RobertUK_F8_clone",   "DMountie",
  "RobertUK_G11_clone",  "DMountie",

  # DMud
  "RobertUK_F3_clone",   "DMud"
)

metadata <- read.csv("/project/berglandlab/Robert/UKSequencing2022_2024/2022_2024seqmetadata20250811.csv", header = TRUE)
# --- apply replacements ---
metadata$Well <- paste0(metadata$Well, "_clone")



metadata4 <- metadata %>%
  mutate(
    accuratelocation_new = if_else(
      is.na(clone) | clone == "" | clone == "BLANK",
      NA_character_,
      str_extract(clone, "^[^_]+")
    )
  )

repl <- metadata4[, c(1, 9)]




metadata2 <- metadata %>%
  left_join(repl, by = "Well") %>%
  mutate(
    accuratelocation = coalesce(accuratelocation_new, accuratelocation)
  ) %>%
  select(-accuratelocation_new)

# sanity check: show Dorset rows after replacement
metadata2 %>%
  filter(grepl("^RobertUK_|^Gilmer5_", clone), grepl("^D", accuratelocation)) %>%
  count(accuratelocation, sort = TRUE)



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
metadata3 <- metadata2 %>%
  filter(!clone %in% c("Blank", "BLANK")) %>%
  mutate(
    clone = paste0(Well),
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
    sampleA = trimws(sample),
    Group   = trimws(mitotype),
    clone   = paste0(sampleA, "_clone")
  ) %>%
  select(clone, mitotype = Group)

## ------------------------
## Tip table: label + mitotype
## Label rule: Sample_ID_old > accuratelocation > clone
## ------------------------



tip_df <- tibble(clone = td@phylo$tip.label) %>%
  left_join(metadata3,  by = "clone") %>%
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
