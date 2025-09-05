
##module load gcc/11.4.0  openmpi/4.1.4 icu R/4.3.1

#R

  library(data.table)
  library(foreach)
  library(SeqArray)
  library(gdsfmt)
library(SNPRelate)
library(ggplot2)
#  library(doMC)
 # registerDoMC(20)

library(ape)
library(dplyr)
library(tibble)
library(igraph)


genofile.fn <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/trimmed10bp_masked_usobtusa.gds"
gds <- seqOpen(genofile.fn)



metadata <- read.csv("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/2022_2024seqmetadata20250811.csv", header = TRUE)
samplestats <- read.csv("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/eudobtusa_samplestats.csv")

metadata_with_clone <- read.csv("/project/berglandlab/Robert/UKSequencing2022_2024/old_stuff/metadata_with_clone.csv", header = TRUE)
metadata_with_clone <- subset(metadata_with_clone, clone !="Blank" & clone !="BLANK")


samples <- read.csv("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/Sample_ID_Species_merged_20250627.csv")
usobtusasamps <- subset(samples, Species == "Daphnia obtusa" & Continent == "NorthAmerica")


samples_to_keep <- usobtusasamps %>% filter(Species == "Daphnia obtusa") %>% pull(Sample_ID)


miss_rate_per_sample <- seqMissing(gds, per.variant = FALSE)
miss_rate_per_variant <- seqMissing(gds, per.variant = TRUE)
sample_ids <- seqGetData(gds, "sample.id")
valid_samples <- sample_ids[miss_rate_per_sample < 0.5]
miss_rate_per_variant <- seqMissing(gds, per.variant=TRUE)
valid_variants <- seqGetData(gds, "variant.id")[miss_rate_per_variant < 0.05]

final_valid_samples <- intersect(valid_samples, samples_to_keep)



seqSetFilter(genofile, sample.id = final_valid_samples)

keep <- which(miss_rate_per_variant < 0.05)



ibs <- snpgdsIBS(gds, sample.id = final_valid_samples, snp.id = keep, num.thread = 10, autosome.only = FALSE)

ibs_mat <- ibs$ibs

# ---- KING IBD (kinship, IBS0, k0/k1 approximation) ----
king <- snpgdsIBDKING(gds, num.thread=10, autosome.only=FALSE)

# Assign row/col names
rownames(king$IBS0) <- king$sample.id
colnames(king$IBS0) <- king$sample.id
rownames(king$kinship) <- king$sample.id
colnames(king$kinship) <- king$sample.id

# Melt both matrices into long form
ibs0_df <- melt(king$IBS0, varnames = c("Sample1", "Sample2"), value.name = "IBS0")
kin_df  <- melt(king$kinship, varnames = c("Sample1", "Sample2"), value.name = "Kinship")

#kin_df <- subset(kin_df, Kinship >= 0)

# Merge
merged_df <- merge(ibs0_df, kin_df, by = c("Sample1", "Sample2"))

# ---- Plot IBS0 vs Kinship ----
p <- ggplot(merged_df, aes(x = IBS0, y = Kinship)) +
  geom_point(alpha = 0.5, size = 1.5, color = "steelblue") +
  theme_bw(base_size = 14) +
  labs(
    title = "IBS0 vs Kinship (KING estimates)",
    x = "IBS0 (proportion sites with 0 shared alleles)",
    y = "Kinship coefficient"
  )

# Save dataframe and plot
write.csv(merged_df, "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/ibs0_kinship_pairs.csv", row.names = FALSE, quote = FALSE)
ggsave("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/ibs0_vs_kinship.png", plot = p, width = 7, height = 6, dpi = 300)


merged_df2 <- left_join(merged_df, metadata, by= c("Sample1"="Well"))



metadata_with_clonekey <- metadata_with_clone[,c(2,6,9)]

merged2 <- left_join(merged_df, metadata_with_clonekey, by = c("Sample1" = "Well"))
colnames(merged2)[5] ="DateA"
colnames(merged2)[6] ="locationA"

merged3 <- left_join(merged2, metadata_with_clonekey, by = c("Sample2" = "Well"))
colnames(merged3)[7] ="DateB"
colnames(merged3)[8] ="locationB"

merged3 <- merged3 %>%
  mutate(
    .a = trimws(as.character(locationA)),
    .b = trimws(as.character(locationB)),
    location_pair = if_else(is.na(.a) | is.na(.b),
                   NA_character_,
                   paste(pmin(.a, .b), pmax(.a, .b), sep = "_"))
  ) %>%
  select(-.a, -.b)

unique(merged3$location_pair)


merged3_sub <- subset(merged3, location_pair == "Gilmer_Gilmer"| 
                      location_pair == "Gilmer_P63"|
                      location_pair == "Gilmer_P66"|
                      location_pair == "Gilmer_P58"|
                      location_pair == "Gilmer_P62"|
                      location_pair == "P63_P63"|
                      location_pair == "P63_P66"|
                      location_pair == "P58_P63"|
                      location_pair == "P62_P63"|
                      location_pair == "P66_P66"|
                      location_pair == "P58_P66"|
                      location_pair == "P62_P66"|
                      location_pair == "P58_P58"|
                      location_pair == "P58_P62"|
                      location_pair == "P62_P62"
    )

# ---- Plot IBS0 vs Kinship ----
p <- ggplot(merged3_sub, aes(x = IBS0, y = Kinship, col=location_pair)) +
  geom_point(alpha = 0.5, size = 1.5) +
  theme_bw(base_size = 14) +
  labs(
    title = "IBS0 vs Kinship (KING estimates)",
    x = "IBS0 (proportion sites with 0 shared alleles)",
    y = "Kinship coefficient"
  )

# Save dataframe and plot
ggsave("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/ibs0_vs_kinship_col_locations.png", plot = p, width = 7, height = 6, dpi = 300)



# ---- Plot IBS0 vs Kinship ----
p <- ggplot(merged3_sub, aes(x = IBS0, y = Kinship, col=location_pair)) +
  geom_point(alpha = 0.5, size = 1.5) +
  theme_bw(base_size = 14) + ylim(0,0.51)+ xlim(0,0.05)+
  labs(
    title = "IBS0 vs Kinship (KING estimates)",
    x = "IBS0 (proportion sites with 0 shared alleles)",
    y = "Kinship coefficient"
  )

# Save dataframe and plot
ggsave("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/ibs0_vs_kinship_col_locations_short.png", plot = p, width = 7, height = 6, dpi = 300)





# ---- Plot IBS0 vs Kinship ----
p <- ggplot(merged3_sub, aes(x = IBS0, y = Kinship, col=location_pair)) +
  geom_point(alpha = 0.5, size = 1.5) +
  theme_bw(base_size = 14) + ylim(0,0.51)+ xlim(0,0.005)+
  labs(
    title = "IBS0 vs Kinship (KING estimates)",
    x = "IBS0 (proportion sites with 0 shared alleles)",
    y = "Kinship coefficient"
  )

# Save dataframe and plot
ggsave("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/ibs0_vs_kinship_col_locations_short_005.png", plot = p, width = 7, height = 6, dpi = 300)






merged4 <- merged3

merged4 <- merged4 %>%
  mutate(
    .a = trimws(as.character(DateA)),
    .b = trimws(as.character(DateB)),
    date_pair = if_else(is.na(.a) | is.na(.b),
                   NA_character_,
                   paste(pmin(.a, .b), pmax(.a, .b), sep = "_"))
  ) %>%
  select(-.a, -.b)

unique(merged4$date_pair)


merged3_sub <- subset(merged3, location_pair == "Gilmer_Gilmer"| 
                      location_pair == "Gilmer_P63"|
                      location_pair == "Gilmer_P66"|
                      location_pair == "Gilmer_P58"|
                      location_pair == "Gilmer_P62"|
                      location_pair == "P63_P63"|
                      location_pair == "P63_P66"|
                      location_pair == "P58_P63"|
                      location_pair == "P62_P63"|
                      location_pair == "P66_P66"|
                      location_pair == "P58_P66"|
                      location_pair == "P62_P66"|
                      location_pair == "P58_P58"|
                      location_pair == "P58_P62"|
                      location_pair == "P62_P62"
    )

date_pair

    "12/22/2023_12/22/2023" 
    "12/22/2023_5/10/2024"  
 "10/7/2023_12/22/2023"                        "12/22/2023_2023"
"12/22/2023_6/8/2021"   "12/22/2023_7/12/2023"  "12/22/2023_5/4/2024"
 "12/22/2023_7/17/2024"  "12/22/2023_6/29/2023"  "12/22/2023_9/2/2023"
 "5/10/2024_5/10/2024"   "2024_5/10/2024"        "10/7/2023_5/10/2024"
 "2023_5/10/2024"        "5/10/2024_6/8/2021"    "5/10/2024_7/12/2023"
 "5/10/2024_5/4/2024"    "5/10/2024_7/17/2024"   "5/10/2024_6/29/2023"
 "5/10/2024_9/2/2023"    "2024_2024"             "10/7/2023_2024"
 "2023_2024"             "2024_6/8/2021"         "2024_7/12/2023"
 "2024_5/4/2024"         "2024_7/17/2024"        "2024_6/29/2023"
 "2024_9/2/2023"         "10/7/2023_10/7/2023"   "10/7/2023_2023"
 "10/7/2023_6/8/2021"    "10/7/2023_7/12/2023"   "10/7/2023_5/4/2024"
 "10/7/2023_7/17/2024"   "10/7/2023_6/29/2023"   "10/7/2023_9/2/2023"
 "2023_2023"             "2023_6/8/2021"         "2023_7/12/2023"
 "2023_5/4/2024"         "2023_7/17/2024"        "2023_6/29/2023"
 "2023_9/2/2023"         "6/8/2021_6/8/2021"     "6/8/2021_7/12/2023"
"5/4/2024_6/8/2021"     "6/8/2021_7/17/2024"    "6/29/2023_6/8/2021"
 "6/8/2021_9/2/2023"     "7/12/2023_7/12/2023"   "5/4/2024_7/12/2023"
 "7/12/2023_7/17/2024"   "6/29/2023_7/12/2023"   "7/12/2023_9/2/2023"
 "5/4/2024_5/4/2024"     "5/4/2024_7/17/2024"    "5/4/2024_6/29/2023"
 "5/4/2024_9/2/2023"     "7/17/2024_7/17/2024"   "6/29/2023_7/17/2024"
 "7/17/2024_9/2/2023"    "6/29/2023_6/29/2023"   "6/29/2023_9/2/2023"
 "9/2/2023_9/2/2023"