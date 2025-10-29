#ijob -A berglandlab -c10 -p standard --mem=40G
#module load gcc/11.4.0  openmpi/4.1.4 icu R/4.3.1
#R

library(data.table)
library(ggplot2)


dxy_raw <- read.table("/scratch/rjp5nc/UK2022_2024/allsites_mito/pixy_mito_whole/pixy_14601_incSRR_dxy.txt",header = TRUE,sep = "\t")


# Clean up
dxy_clean <- na.omit(dxy_raw)

dxy_clean$avg_dxy   <- as.numeric(dxy_clean$avg_dxy)

mu_min <- 1.37e-7   # lower rate
mu_max <- 1.73e-7   # upper rate

# --- Calculate divergence times (haploid, so T = dxy / (2 * mu)) ---
dxy_clean$T_min_gen <- dxy_clean$avg_dxy / (2 * mu_max)   # generations, younger bound
dxy_clean$T_max_gen <- dxy_clean$avg_dxy / (2 * mu_min)   # generations, older bound

# --- Convert to years (6 generations per year) ---
dxy_clean$T_min_yrs <- dxy_clean$T_min_gen / 6
dxy_clean$T_max_yrs <- dxy_clean$T_max_gen / 6


# Save plot as PNG
png("/scratch/rjp5nc/UK2022_2024/allsites_mito/dxy_hist.png",
    width = 3000, height = 2000, res = 300)
ggplot(dxy_clean, aes(x = avg_dxy)) +
  geom_histogram(color = "black", fill = "steelblue", bins = 30) +
  theme_classic(base_size = 14) +
  labs(
    x = expression("Average D"[XY]),
    y = "Number of pairwise comparisons",
    title = "Distribution of average pairwise divergence (dxy)"
  ) 
dev.off()



genomic_types <- read.csv("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/genomic_types.csv")
mitotypes <- read.csv("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/mito_types.csv")


library(dplyr)

# Rename columns
genomic_types <- genomic_types %>%
  rename(superclone = Group,
         sample = CloneA)

mitotypes <- mitotypes %>%
  rename(mitotype = Group,
         sample = CloneA)

# Join both annotations to dxy_clean
dxy_annot <- dxy_clean %>%
  # add superclone info for pop1 and pop2
  left_join(genomic_types, by = c("pop1" = "sample")) %>%
  rename(superclone_pop1 = superclone) %>%
  left_join(genomic_types, by = c("pop2" = "sample")) %>%
  rename(superclone_pop2 = superclone) %>%
  # add mitotype info for pop1 and pop2
  left_join(mitotypes, by = c("pop1" = "sample")) %>%
  rename(mitotype_pop1 = mitotype) %>%
  left_join(mitotypes, by = c("pop2" = "sample")) %>%
  rename(mitotype_pop2 = mitotype)

# Check the first few rows
head(dxy_annot)

dxy_filtered <- dxy_annot %>%
  filter(
    # Keep if both annotations exist
    (!is.na(superclone_pop1) & !is.na(superclone_pop2) &
     !is.na(mitotype_pop1) & !is.na(mitotype_pop2)) |
    # OR keep if either pop1 or pop2 starts with "SRR"
    grepl("^SRR", pop1) | grepl("^SRR", pop2)
  )

  dxy_filtered <- dxy_filtered %>%
  mutate(
    superclone_pop1 = ifelse(grepl("^SRR", pop1), "XX", superclone_pop1),
    superclone_pop2 = ifelse(grepl("^SRR", pop2), "XX", superclone_pop2),
    mitotype_pop1   = ifelse(grepl("^SRR", pop1), "XX", mitotype_pop1),
    mitotype_pop2   = ifelse(grepl("^SRR", pop2), "XX", mitotype_pop2)
  )





dxy_filtered$mito_comparison <- paste(dxy_filtered$mitotype_pop1, dxy_filtered$mitotype_pop2, sep = "_")
unique(dxy_filtered$mito_comparison)

dxy_filtered <- dxy_filtered %>%
  mutate(
    mito_comparison = paste(mitotype_pop1, mitotype_pop2, sep = "_"),
    mito_pair = apply(
      cbind(mitotype_pop1, mitotype_pop2),
      1,
      function(x) paste(sort(x), collapse = "_")
    )
  )


missing_mito <- dxy_filtered %>%
  filter(is.na(mitotype_pop1) | is.na(mitotype_pop2) |
         mitotype_pop1 == "XX" | mitotype_pop2 == "XX")


dxy_filtered <- dxy_filtered %>%
  filter(
    # Keep if both annotations exist
    (!is.na(superclone_pop1) & !is.na(superclone_pop2) &
     !is.na(mitotype_pop1) & !is.na(mitotype_pop2)))


dxy_filtered <- subset(dxy_filtered, count_comparisons >= 500)

png("/scratch/rjp5nc/UK2022_2024/allsites_mito/dxy_hist_self.png",
    width = 3000, height = 2000, res = 300)
ggplot(subset(dxy_filtered, mito_pair %in% c("A_A", "B_B", "C_C", "D_D", "E_E")), aes(x = avg_dxy)) +
  geom_histogram(color = "black", fill = "steelblue", bins = 30) +
  theme_classic(base_size = 14) +
    facet_wrap(~ mito_pair, scales = "free_y") +
  labs(
    x = expression("Average D"[XY]),
    y = "Number of pairwise comparisons",
    title = "Distribution of average pairwise divergence (dxy)"
  ) 
dev.off()

png("/scratch/rjp5nc/UK2022_2024/allsites_mito/dxy_hist_nonself_no_ef.png",
    width = 3000, height = 2000, res = 300)

ggplot(
  subset(
    dxy_filtered,
    !(mito_pair %in% c("A_A", "B_B", "C_C", "D_D", "E_E")) & !grepl("E", mito_pair) & !grepl("F", mito_pair)
  ),
  aes(x = avg_dxy)
) +
  geom_histogram(color = "black", fill = "steelblue", bins = 30) +
  theme_classic(base_size = 14) +
  facet_wrap(~ mito_pair, scales = "free_y") +
  labs(
    x = expression("Average D"[XY]),
    y = "Number of pairwise comparisons",
    title = "Distribution of average pairwise divergence (dxy) between mitotypes"
  )

dev.off()











dxy_filtered <- dxy_filtered %>%
  mutate(
    super_comparison = paste(superclone_pop1, superclone_pop2, sep = "_"),
    super_pair = apply(
      cbind(superclone_pop1, superclone_pop2),
      1,
      function(x) paste(sort(x), collapse = "_")
    )
  )




  
png("/scratch/rjp5nc/UK2022_2024/allsites_mito/dxy_hist_super_self.png",
    width = 5000, height = 5000, res = 300)
ggplot(subset(dxy_filtered, superclone_pop1 == superclone_pop2), aes(x = avg_dxy)) +
  geom_histogram(color = "black", fill = "steelblue", bins = 30) +
  theme_classic(base_size = 14) +
    facet_wrap(~ super_pair, scales = "free_y") +
  labs(
    x = expression("Average D"[XY]),
    y = "Number of pairwise comparisons",
    title = "Distribution of average pairwise divergence (dxy)"
  ) 
dev.off()


metadata_with_clone <- read.csv("/project/berglandlab/Robert/UKSequencing2022_2024/old_stuff/metadata_with_clone.csv", header = TRUE)
metadata_with_clone$clone <- trimws(metadata_with_clone$clone)
metadata_with_clone <- subset(metadata_with_clone, clone !="Blank")
metadata_with_clone <- subset(metadata_with_clone, clone !="BLANK")



metadata_subset <- metadata_with_clone[, c("Well", "clone", "accuratelocation")]
head(metadata_subset)


dxy_loc <- dxy_filtered %>%
  # Join location info for pop1
  left_join(metadata_subset, by = c("pop1" = "Well")) %>%
  rename(clone_pop1 = clone,
         location_pop1 = accuratelocation) %>%
  # Join location info for pop2
  left_join(metadata_subset, by = c("pop2" = "Well")) %>%
  rename(clone_pop2 = clone,
         location_pop2 = accuratelocation) %>%
  # Assign VermilionFairgrounds to SRR samples
  mutate(
    location_pop1 = ifelse(grepl("^SRR", pop1), "VermilionFairgrounds", location_pop1),
    location_pop2 = ifelse(grepl("^SRR", pop2), "VermilionFairgrounds", location_pop2)
  )



dxy_loc2 <- dxy_loc %>%
  mutate(
    location_comparison = paste(location_pop1, location_pop2, sep = "_"),
    location_pair = apply(
      cbind(location_pop1, location_pop2),
      1,
      function(x) paste(sort(x), collapse = "_")
    )
  )





png("/scratch/rjp5nc/UK2022_2024/allsites_mito/dxy_hist_loc_self.png",
    width = 3000, height = 2000, res = 300)
ggplot(subset(dxy_loc2, location_pop1 == location_pop2), aes(x = avg_dxy)) +
  geom_histogram(color = "black", fill = "steelblue", bins = 30) +
  theme_classic(base_size = 14) +
    facet_wrap(~ location_pair, scales = "free_y") +
  labs(
    x = expression("Average D"[XY]),
    y = "Number of pairwise comparisons",
    title = "Distribution of average pairwise divergence (dxy)"
  ) 
dev.off()


png("/scratch/rjp5nc/UK2022_2024/allsites_mito/dxy_hist_loc_nonself.png",
    width = 3000, height = 2000, res = 300)

ggplot(
  subset(dxy_loc2, location_pop1 != location_pop2),
  aes(x = avg_dxy)
) +
  geom_histogram(color = "black", fill = "steelblue", bins = 30) +
  theme_classic(base_size = 14) +
  facet_wrap(~ location_pair, scales = "free_y") +
  labs(
    x = expression("Average D"[XY]),
    y = "Number of pairwise comparisons",
    title = "Distribution of average pairwise divergence (dxy) between mitotypes"
  )

dev.off()






div_times <- dxy_loc2 %>%
  mutate(
    T_mid_yrs = (T_min_yrs + T_max_yrs) / 2,
    T_err_yrs = (T_max_yrs - T_min_yrs) / 2
  )

png("/scratch/rjp5nc/UK2022_2024/allsites_mito/dxy_hist_div.png",
    width = 3000, height = 2000, res = 300)

ggplot(
  subset(div_times, !is.na(T_mid_yrs) & !is.na(mito_pair)),
  aes(x = mito_pair, y = T_mid_yrs)
) +
  geom_jitter(width = 0.15, height = 0, alpha = 0.35) +
  geom_boxplot(fill = NA, color = "black", outlier.shape = NA) +
  theme_classic(base_size = 14) +
  labs(
    x = "Mitotype pair",
    y = "Divergence time (years)",
    title = "Divergence times between lineages by mitotype pair"
  )
  dev.off()






png("/scratch/rjp5nc/UK2022_2024/allsites_mito/dxy_hist_div_noSRR.png",
    width = 3000, height = 2000, res = 300)

ggplot(
  subset(div_times, !is.na(T_mid_yrs) & !is.na(mito_pair) & !grepl("XX", mito_pair)),
  aes(x = mito_pair, y = T_mid_yrs)
) +
  geom_jitter(width = 0.15, height = 0, alpha = 0.35) +
  geom_boxplot(fill = NA, color = "black", outlier.shape = NA) +
  theme_classic(base_size = 14) +
  labs(
    x = "Mitotype pair",
    y = "Divergence time (years)",
    title = "Divergence times between lineages by mitotype pair (no SRR/XX)"
  )

dev.off()





png("/scratch/rjp5nc/UK2022_2024/allsites_mito/dxy_hist_div_major3.png",
    width = 3000, height = 2000, res = 300)

ggplot(
  subset(div_times, !is.na(T_mid_yrs) & !is.na(mito_pair) & !grepl("XX", mito_pair) & !grepl("E", mito_pair) & !grepl("F", mito_pair)),
  aes(x = mito_pair, y = T_mid_yrs, color=location_pair)
) +
  geom_jitter(width = 0.15, height = 0, alpha = 0.35) +
  geom_boxplot(fill = NA, color = "black", outlier.shape = NA) +
  theme_classic(base_size = 14) +
  labs(
    x = "Mitotype pair",
    y = "Divergence time (years)",
    title = "Divergence times between lineages by mitotype pair (no SRR/XX)"
  )

dev.off()





png("/scratch/rjp5nc/UK2022_2024/allsites_mito/dxy_hist_div_major3_SRR.png",
    width = 5000, height = 2000, res = 300)

ggplot(
  subset(div_times, !is.na(T_mid_yrs) & !is.na(mito_pair) & !grepl("E", mito_pair) & !grepl("F", mito_pair)),
  aes(x = mito_pair, y = T_mid_yrs, color=location_pair)
) +
  geom_jitter(width = 0.15, height = 0, alpha = 0.35) +
  geom_boxplot(fill = NA, color = "black", outlier.shape = NA) +
  theme_classic(base_size = 14) +
  labs(
    x = "Mitotype pair",
    y = "Divergence time (years)",
    title = "Divergence times between lineages by mitotype pair (no SRR/XX)"
  )

dev.off()


subset(div_times, mito_pair == "A_XX" & T_mid_yrs <= 10000)

summary(div_times$count_comparisons)


# # Save plot as PNG
# png("/scratch/rjp5nc/UK2022_2024/allsites_mito/counts.png",
#     width = 3000, height = 2000, res = 300)
# ggplot(div_times, aes(x = count_comparisons)) +
#   geom_histogram(color = "black", fill = "steelblue", bins = 30) +
#   theme_classic(base_size = 14) +
#   labs(
#     x = expression("Average D"[XY]),
#     y = "Count comparisons",
#     title = "count comparisons"
#   ) 
# dev.off()




summary(subset(div_times, count_diffs == 1)$T_mid_yrs)





png("/scratch/rjp5nc/UK2022_2024/allsites_mito/dxy_hist_div_AA.png",
    width = 5000, height = 2000, res = 300)

ggplot(
  subset(div_times, mito_pair == "A_A"),
  aes(x = location_pair, y = T_mid_yrs, color=location_pair)
) +
  geom_jitter(width = 0.15, height = 0, alpha = 0.35) +
  geom_boxplot(fill = NA, color = "black", outlier.shape = NA) +
  theme_classic(base_size = 14) +
  labs(
    x = "Mitotype pair",
    y = "Divergence time (years)",
    title = "Divergence times between lineages by mitotype pair (A_A)"
  )

dev.off()
