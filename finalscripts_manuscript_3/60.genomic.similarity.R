##module load gcc/11.4.0  openmpi/4.1.4 icu R/4.3.1

#R

# Load required packages
library(SeqArray)
library(SNPRelate)
library(ape)
library(ggtree)
library(ggplot2)
library(dplyr)
library(tibble)
library(data.table)
library(igraph)

similarity <- read.csv("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/usobtusa_similarity_matrix.csv")


metadata_with_clone <- read.csv("/project/berglandlab/Robert/UKSequencing2022_2024/old_stuff/metadata_with_clone.csv", header = TRUE)
samples <- read.csv("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/Sample_ID_Species_merged_20250627.csv")
usobtusasamps <- subset(samples, Species == "Daphnia obtusa" & Continent == "NorthAmerica")
metadata_with_clone <- subset(metadata_with_clone, clone !="Blank" & clone !="BLANK")

# turn rownames back into row names
rownames(similarity) <- similarity$X
similarity$X <- NULL

# now convert to matrix and melt
similarity_mat <- as.matrix(similarity)

long_mat <- reshape2::melt(similarity_mat, varnames = c("CloneA", "CloneB"), value.name = "Similarity") %>%
  mutate(CloneA = as.character(CloneA),
         CloneB = as.character(CloneB))

long_mat$Similarity <- 1 - long_mat$Similarity


metadata_with_clonekey <- metadata_with_clone[,c(2,9)]

long_matmerged1 <- left_join(long_mat, metadata_with_clonekey, by = c("CloneA" = "Well"))
colnames(long_matmerged1)[4] ="locationA"

long_matmerged2 <- left_join(long_matmerged1, metadata_with_clonekey, by = c("CloneB" = "Well"))
colnames(long_matmerged2)[5] ="locationB"

long_matmerged2 <- long_matmerged2 %>%
  mutate(
    .a = trimws(as.character(locationA)),
    .b = trimws(as.character(locationB)),
    pair = if_else(is.na(.a) | is.na(.b),
                   NA_character_,
                   paste(pmin(.a, .b), pmax(.a, .b), sep = "_"))
  ) %>%
  select(-.a, -.b)

pairs_sorted <- sort(unique(na.omit(long_matmerged2$pair)))

# create a mapping of pair -> letter
pair_letters <- setNames(LETTERS[seq_along(pairs_sorted)], pairs_sorted)

# assign letters
long_matmerged2 <- long_matmerged2 %>%
  mutate(group_letter = pair_letters[pair])

unique(long_matmerged2$group_letter)
unique(long_matmerged2$pair)

grouplsplots <- ggplot(long_matmerged2, aes(x = CloneA, y=CloneB, col=Similarity)) +
geom_point()+
  theme_bw() + 
  xlab("Similarity") +
  ylab("Frequency") 

ggsave("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/groupsplot.png", plot = grouplsplots, width = 12, height = 8, dpi = 300)

pair_means <- long_matmerged2 %>%
  group_by(pair) %>%
  summarise(mean_similarity = mean(Similarity, na.rm = TRUE)) %>%
  arrange(pair)

as.data.frame(pair_means)

pair_stats <- long_matmerged2 %>%
  group_by(pair) %>%
  summarise(
    mean_similarity = mean(Similarity, na.rm = TRUE),
    median_similarity = median(Similarity, na.rm = TRUE)
  ) %>%
  arrange(pair)

as.data.frame(pair_stats)