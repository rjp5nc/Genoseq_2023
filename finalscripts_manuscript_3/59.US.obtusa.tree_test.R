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

# ---- Step 1: Open GDS file ----

metadata <- read.csv("/project/berglandlab/Robert/UKSequencing2022_2024/old_stuff/2022_2024seqmetadata20250131.csv", header = TRUE)
metadata_with_clone <- read.csv("/project/berglandlab/Robert/UKSequencing2022_2024/old_stuff/metadata_with_clone.csv", header = TRUE)
samples <- read.csv("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/Sample_ID_Species_merged_20250627.csv")
usobtusasamps <- subset(samples, Species == "Daphnia obtusa" & Continent == "NorthAmerica")
metadata_with_clone <- subset(metadata_with_clone, clone !="Blank" & clone !="BLANK")
mitotypes <- read.csv("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/mito_types.csv")


head(metadata)
seqClose(genofile)  
genofile.fn <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/trimmed10bp_masked_usobtusa.gds"
genofile <- seqOpen(genofile.fn)


seqResetFilter(genofile)

samples_to_keep <- usobtusasamps %>% filter(Species == "Daphnia obtusa") %>% pull(Sample_ID)
#samples_to_keep <- metadata_with_clone %>% filter(accuratelocation == "P58") %>% pull(Well)
#samples_to_keep <- metadata_with_clone %>% filter(location == "UK" | accuratelocation == "P759") %>% pull(Well)

#samples_to_keep <- metadata_with_clone$Well

unique(metadata_with_clone$accuratelocation)

# ---- Step 2: Filter variants with missing rate < 0.05 ----
miss_rate_per_sample <- seqMissing(genofile, per.variant = FALSE)
miss_rate_per_variant <- seqMissing(genofile, per.variant = TRUE)
sample_ids <- seqGetData(genofile, "sample.id")
valid_samples <- sample_ids[miss_rate_per_sample < 0.5]
miss_rate_per_variant <- seqMissing(genofile, per.variant=TRUE)
valid_variants <- seqGetData(genofile, "variant.id")[miss_rate_per_variant < 0.05]

final_valid_samples <- intersect(valid_samples, samples_to_keep)


seqSetFilter(genofile, sample.id = final_valid_samples)

miss_rate <- seqMissing(genofile, per.variant = TRUE)
dp <- seqGetData(genofile, "annotation/format/DP")
mean_depth <- rowMeans(dp, na.rm = TRUE)
keep <- which(miss_rate < 0.05)


# ---- Step 3: Compute IBS distance matrix ----
ibs <- snpgdsIBS(genofile, sample.id = final_valid_samples, snp.id = keep, num.thread = 4, autosome.only = FALSE)
dist_matrix <- 1 - ibs$ibs
rownames(dist_matrix) <- colnames(dist_matrix) <- ibs$sample.id
write.csv(dist_matrix, "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/usobtusa_similarity_matrix.csv")  # Save as CSV

dist_matrix2 <- 1- dist_matrix

long_mat <- melt(dist_matrix2, varnames = c("CloneA", "CloneB"), value.name = "Similarity") %>%
  mutate(CloneA = as.character(CloneA),
         CloneB = as.character(CloneB))  # ensure they are character, not factor


long_mat_plot <- ggplot(long_mat, aes(x = Similarity)) +
  geom_histogram(bins = 200) +
  theme_bw() + 
  xlab("Similarity") + ylim(0,15000) + xlim(0.5,1.05) + 
  ylab("Frequency") +
  geom_vline(xintercept = 0.975, linetype = "dashed", color = "red") 

ggsave("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/long_mat_genomic_plot.png", plot = long_mat_plot, width = 6, height = 4, dpi = 300)

# 2. Filter dataset
long_filt <- long_mat %>%
  filter(Similarity >= 0.965, Similarity != 1)

# 3. Build graph
g <- graph_from_data_frame(long_filt[, c("CloneA", "CloneB")], directed = FALSE)


png("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/grouping_graph_genomic.png",
    width = 2000, height = 2000, res = 300)
plot(g, vertex.size = 5, vertex.label = NA)  # basic igraph plot
dev.off()

# 4. Find connected components
comp <- components(g)


labels <- c(LETTERS, paste0("A", LETTERS), paste0("B", LETTERS))

# Map them to your groups
group_letters <- setNames(labels[comp$membership], names(comp$membership))

# 5. Assign group letters

# 6. Add group column (both clones exist in the graph)
long_filt <- long_filt %>%
  mutate(Group = group_letters[CloneA])


head(long_filt)
unique(long_filt$Group)

long_filt_onlywell <- long_filt[,c(1,4)]

unique_clones <- long_filt_onlywell %>%
  distinct(CloneA, Group)
unique_clones
unique(long_filt$Group)


# ---- Step 4: Build neighbor-joining tree ----
tree <- nj(as.dist(dist_matrix))

# ---- Step 5 (Optional): Relabel tips from metadata CSV ----
# CSV should have columns: Sample_ID and Label (or Well, Clone, etc.)
tree_rooted <- root(tree, outgroup = "SRR5012393", resolve.root = TRUE)

merged_data2 <- left_join(tree_rooted, metadata_with_clone, by = c("label" = "Well"))


group_colors <- setNames(rainbow(length(unique(unique_clones$Group))), unique(unique_clones$Group))

# Assign color to each tip

tip_colors <- group_colors[unique_clones$Group]

tip_color_dt <- data.table(unique_clones$Group, unique_clones$CloneA)

final_valid_samples2 <- as.data.frame(final_valid_samples)

final_valid_samples3 <- left_join(final_valid_samples2, tip_color_dt, by= c("final_valid_samples"="V2"))

tip_colors <- group_colors[final_valid_samples3$V1]


tip_color_dt <- data.table(
  sample = names(tip_colors),
  color  = as.vector(tip_colors)
)

phylo_tree <- merged_data2@phylo
write.tree(phylo_tree, file = "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/usdobtusa_tree_genomic.nwk")

png("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/usobtusa_groups.png", res = 300, width = 4000, height = 9000)
#png("/scratch/rjp5nc/UK2022_2024/mito_vcf/tree_usobtusa_circ.png", width = 1200, height = 2000)

# Plot the tree
plot.phylo(phylo_tree,
           type = "phylogram",
#           type = "fan",
           cex = 0.8,
           label.offset = 0.01,
           no.margin = TRUE,
           tip.color = tip_colors,
           main = "Rooted NJ Tree (root = P759)")

# Add legend
legend("topleft",                   # position
       legend = unique(tip_color_dt$sample),     # group names
       col = unique(tip_color_dt$color),        # matching colors
       pch = 19,                   # solid circle
       pt.cex = 1.5,               # point size
       cex = 1,                    # text size
       bty = "n",                  # no box
       title = "Sample Group")

dev.off()

write.csv(unique_clones, "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/genomic_types.csv")  # Save as CSV





merged_data3 <- left_join(merged_data2, mitotypes, by = c("label" = "CloneA"))


p <- ggtree(merged_data3, options(ignore.negative.edge = TRUE)) +
  geom_tiplab(aes(label = clone, color = Group), size = 3) +   # <- use clone + location
  theme_tree() +
  labs(title = "Rooted NJ Tree root SRR") 

ggsave(
  "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/usobtusa_groups2.png",
  plot = p,
  width = 10, height = 40, dpi = 300
)



p <- ggtree(merged_data3, options(ignore.negative.edge = TRUE), layout="circular") +
  geom_tiplab(aes(label = clone, color = Group), size = 3) +   # <- use clone + location
  theme_tree() +
  labs(title = "Rooted NJ Tree root SRR") 

ggsave(
  "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/usobtusa_groups2_fan.png",
  plot = p,
  width = 20, height = 20, dpi = 300
)





# ---- Step 4: Build neighbor-joining tree ----
tree <- nj(as.dist(dist_matrix))

# ---- Step 5 (Optional): Relabel tips from metadata CSV ----
# CSV should have columns: Sample_ID and Label (or Well, Clone, etc.)
tree_rooted <- root(tree, outgroup = "SRR5012394", resolve.root = TRUE)

label_map <- setNames(mitotypes$CloneA, mitotypes$Group)

# Apply labels (optional: fallback to Sample_ID if label is missing)
tree$tip.label <- label_map[tree$tip.label]
tree$tip.label[is.na(tree$tip.label)] <- names(label_map)[is.na(tree$tip.label)]

group_colors <- setNames(rainbow(length(unique(mitotypes$Group))), unique(mitotypes$Group))

# Assign color to each tip

tip_colors <- group_colors[mitotypes$Group]

tip_color_dt <- data.table(mitotypes$Group, mitotypes$CloneA)

final_valid_samples2 <- as.data.frame(final_valid_samples)

final_valid_samples3 <- left_join(final_valid_samples2, tip_color_dt, by= c("final_valid_samples"="V2"))

tip_colors <- group_colors[final_valid_samples3$V1]


tip_color_dt <- data.table(
  sample = names(tip_colors),
  color  = as.vector(tip_colors)
)


png("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/usobtusa_genomicxmitogroups.png", res = 300, width = 4000, height = 9000)
#png("/scratch/rjp5nc/UK2022_2024/mito_vcf/tree_usobtusa_circ.png", width = 1200, height = 2000)

# Plot the tree
plot.phylo(tree_rooted,
           type = "phylogram",
#           type = "fan",
           cex = 0.8,
           label.offset = 0.01,
           no.margin = TRUE,
           tip.color = tip_colors,
           main = "Rooted NJ Tree (root = P759)")

# Add legend
legend("topleft",                   # position
       legend = unique(tip_color_dt$sample),     # group names
       col = unique(tip_color_dt$color),        # matching colors
       pch = 19,                   # solid circle
       pt.cex = 1.5,               # point size
       cex = 1,                    # text size
       bty = "n",                  # no box
       title = "Sample Group")

dev.off()









# ---- Step 4: Build neighbor-joining tree ----
tree <- nj(as.dist(dist_matrix))

# ---- Step 5 (Optional): Relabel tips from metadata CSV ----
# CSV should have columns: Sample_ID and Label (or Well, Clone, etc.)
tree_rooted <- root(tree, outgroup = "Gilmer5_E2", resolve.root = TRUE)

label_map <- setNames(metadata_with_clone$clone, metadata_with_clone$Well)

# Apply labels (optional: fallback to Sample_ID if label is missing)
tree$tip.label <- label_map[tree$tip.label]
tree$tip.label[is.na(tree$tip.label)] <- names(label_map)[is.na(tree$tip.label)]

group_colors <- setNames(rainbow(length(unique(metadata_with_clone$date))), unique(metadata_with_clone$date))

# Assign color to each tip

tip_colors <- group_colors[metadata_with_clone$date]

tip_color_dt <- data.table(metadata_with_clone$date, metadata_with_clone$Well)

final_valid_samples2 <- as.data.frame(final_valid_samples)

final_valid_samples3 <- left_join(final_valid_samples2, tip_color_dt, by= c("final_valid_samples"="V2"))

tip_colors <- group_colors[final_valid_samples3$V1]


tip_color_dt <- data.table(
  sample = names(tip_colors),
  color  = as.vector(tip_colors)
)


png("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/usobtusa_genomic.png", res = 300, width = 4000, height = 9000)
#png("/scratch/rjp5nc/UK2022_2024/mito_vcf/tree_usobtusa_circ.png", width = 1200, height = 2000)

# Plot the tree
plot.phylo(tree_rooted,
           type = "phylogram",
#           type = "fan",
           cex = 0.8,
           label.offset = 0.01,
           no.margin = TRUE,
           tip.color = tip_colors,
           main = "Rooted NJ Tree (root = P759)")

# Add legend
legend("topleft",                   # position
       legend = unique(tip_color_dt$sample),     # group names
       col = unique(tip_color_dt$color),        # matching colors
       pch = 19,                   # solid circle
       pt.cex = 1.5,               # point size
       cex = 1,                    # text size
       bty = "n",                  # no box
       title = "Sample Group")

dev.off()




png("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/usobtusa_genomic_rooted_by_SRR.png", res = 300, width = 4000, height = 4000)
#png("/scratch/rjp5nc/UK2022_2024/mito_vcf/tree_usobtusa_circ.png", width = 1200, height = 2000)

# Plot the tree
plot.phylo(tree_rooted,
           type = "fan",
           cex = 0.8,
           label.offset = 0.01,
           no.margin = TRUE,
           tip.color = tip_colors,
           main = "Rooted NJ Tree")

# Add legend
legend("topleft",                   # position
       legend = unique(tip_color_dt$sample),     # group names
       col = unique(tip_color_dt$color),        # matching colors
       pch = 19,                   # solid circle
       pt.cex = 1.5,               # point size
       cex = 1,                    # text size
       bty = "n",                  # no box
       title = "Sample Group")

dev.off()


metadata_sub <- metadata %>%
  filter(Well %in% final_valid_samples)
metadata_sub <- metadata_sub[metadata_sub$accuratelocation != "PBO66", ]

#metadata_sub <- subset(metadata_sub, accuratelocation == "P62"| accuratelocation == "P66" )

unique(metadata_sub$accuratelocation)

pop_factor <- factor(metadata_sub$accuratelocation)

seqResetFilter(genofile)

seqSetFilter(genofile, sample.id = metadata_sub$Well)

write.csv(metadata_sub, "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/metadata_sub.csv")  # Save as CSV

metadata_sub_sum <- metadata_sub %>%
  group_by(accuratelocation, date) %>%
  summarise(count = n(), .groups = "drop")

fst <- snpgdsFst(
  genofile,
  sample.id = metadata_sub$Well,
  population = pop_factor,
  snp.id = valid_variants,
  autosome.only = FALSE,
  method = "W&C84",
  maf= 0.05,
  missing.rate= 0.5, verbose=TRUE
)

fst$Fst





#FST overall between ponds

pop_levels <- unique(metadata_sub$accuratelocation)
pairwise_results <- matrix(NA, nrow=length(pop_levels), ncol=length(pop_levels),
                           dimnames=list(pop_levels, pop_levels))

for(i in 1:(length(pop_levels)-1)){
  for(j in (i+1):length(pop_levels)){
    # logical index for samples in populations i and j
    idx <- metadata_sub$accuratelocation %in% c(pop_levels[i], pop_levels[j])
    
    sub_samples <- metadata_sub$Well[idx]
    sub_pop <- factor(metadata_sub$accuratelocation[idx])
    
    # Only calculate if both populations exist
    if(length(unique(sub_pop)) == 2){
      fst_tmp <- snpgdsFst(
        genofile,
        sample.id = sub_samples,
        population = sub_pop,
        snp.id = valid_variants,
        autosome.only = FALSE,
        method = "W&C84",
        maf = 0.05,
        missing.rate = 0.5,
        verbose = FALSE
      )
      
      pairwise_results[i,j] <- fst_tmp$Fst
      pairwise_results[j,i] <- fst_tmp$Fst
    }
  }
}

# Fill diagonal
diag(pairwise_results) <- 0
pairwise_results

write.csv(pairwise_results,"/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/pairwise_fst.csv")

# Melt the matrix
fst_long <- melt(pairwise_results, varnames = c("PopulationA", "PopulationB"), value.name = "Fst")

# Convert to character to compare
fst_long$PopulationA <- as.character(fst_long$PopulationA)
fst_long$PopulationB <- as.character(fst_long$PopulationB)

# Remove duplicates and diagonal
fst_long <- fst_long[fst_long$PopulationA < fst_long$PopulationB, ]

fst_long


fst_melted <- melt(pairwise_results, varnames = c("PopulationA", "PopulationB"), value.name = "Fst")

# Heatmap
popheatmapplot <- ggplot(fst_melted, aes(x = PopulationA, y = PopulationB, fill = Fst)) +
  geom_tile(color = "white") +
  geom_text(aes(label = ifelse(is.na(Fst), "", sprintf("%.3f", Fst))), size = 3) +
  scale_fill_gradient(low = "white", high = "red", na.value = "grey90") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    panel.grid = element_blank()
  ) +
  labs(title = "Pairwise Fst Heatmap", fill = "Fst")


ggsave("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/popheatmap.png", plot = popheatmapplot, width = 7, height = 6, dpi = 300)












#FST by each date

library(dplyr)


metadata_sub <- metadata %>%
  filter(Well %in% final_valid_samples)
metadata_sub <- metadata_sub[metadata_sub$accuratelocation != "PBO66", ]

#metadata_sub <- subset(metadata_sub, accuratelocation == "P62"| accuratelocation == "P66" )

unique(metadata_sub$accuratelocation)

seqResetFilter(genofile)

seqSetFilter(genofile, sample.id = metadata_sub$Well)

# unique dates in your metadata
date_levels <- unique(metadata_sub$date)
pop_levels  <- unique(metadata_sub$accuratelocation)

# create a list to hold fst matrices per date
fst_by_date <- list()

for (d in date_levels) {
  # subset metadata for this date
  meta_d <- metadata_sub %>% filter(date == d)
  
  # initialize empty matrix for this date
  pairwise_results <- matrix(NA, 
                             nrow = length(pop_levels), 
                             ncol = length(pop_levels),
                             dimnames = list(pop_levels, pop_levels))
  
  for (i in 1:(length(pop_levels) - 1)) {
    for (j in (i + 1):length(pop_levels)) {
      
      idx <- meta_d$accuratelocation %in% c(pop_levels[i], pop_levels[j])
      sub_samples <- meta_d$Well[idx]
      sub_pop <- factor(meta_d$accuratelocation[idx])
      
      if (length(unique(sub_pop)) == 2) {
        fst_tmp <- snpgdsFst(
          genofile,
          sample.id = sub_samples,
          population = sub_pop,
          snp.id = valid_variants,
          autosome.only = FALSE,
          method = "W&C84",
          maf = 0.05,
          missing.rate = 0.5,
          verbose = FALSE
        )
        
        pairwise_results[i, j] <- fst_tmp$Fst
        pairwise_results[j, i] <- fst_tmp$Fst
      }
    }
  }
  
  fst_by_date[[d]] <- pairwise_results
}

library(dplyr)

# flatten the list into a long dataframe
fst_long_date <- do.call(rbind, lapply(names(fst_by_date), function(d) {
  mat <- fst_by_date[[d]]
  if (is.null(mat)) return(NULL)
  
  df <- as.data.frame(as.table(mat)) %>%
    rename(pop1 = Var1, pop2 = Var2, Fst = Freq) %>%
    mutate(date = d)
  
  return(df)
}))

# ensure pop1/pop2 are characters
fst_long_date <- fst_long_date %>%
  mutate(pop1 = as.character(pop1),
         pop2 = as.character(pop2))

# remove diagonal and duplicates
fst_long_date <- fst_long_date %>%
  filter(!is.na(Fst), pop1 < pop2)

# write to CSV
write.csv(fst_long_date,
          "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/fst_by_date.csv",
          row.names = FALSE)


bydateplot <- ggplot(fst_long_date, aes(x = pop1, y = pop2, fill = Fst)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.3f", Fst)), size = 3) +
  scale_fill_gradient(low = "white", high = "red", na.value = "grey90") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    panel.grid = element_blank()
  ) +
  labs(title = "Pairwise Fst by Date", fill = "Fst") +
  facet_wrap(~ date)


ggsave("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/popheatmap_bydateplot.png", plot = bydateplot, width = 7, height = 6, dpi = 300)






#FST within the same pond by date

# unique ponds
ponds <- unique(metadata_sub$accuratelocation)

# list to hold results
fst_within_pond <- list()

for (pond in ponds) {
  meta_p <- metadata_sub %>% filter(accuratelocation == pond)
  date_levels <- unique(meta_p$date)
  
  # matrix of pairwise Fst between dates for this pond
  pairwise_results <- matrix(NA,
                             nrow = length(date_levels),
                             ncol = length(date_levels),
                             dimnames = list(date_levels, date_levels))
  
  for (i in 1:(length(date_levels) - 1)) {
    for (j in (i + 1):length(date_levels)) {
      
      idx <- meta_p$date %in% c(date_levels[i], date_levels[j])
      sub_samples <- meta_p$Well[idx]
      sub_pop <- factor(meta_p$date[idx])  # population defined by date
      
      # only compute if both dates have samples
      if (length(unique(sub_pop)) == 2) {
        fst_tmp <- snpgdsFst(
          genofile,
          sample.id = sub_samples,
          population = sub_pop,
          snp.id = valid_variants,
          autosome.only = FALSE,
          method = "W&C84",
          maf = 0.05,
          missing.rate = 0.5,
          verbose = FALSE
        )
        
        pairwise_results[i, j] <- fst_tmp$Fst
        pairwise_results[j, i] <- fst_tmp$Fst
      }
    }
  }
  
  fst_within_pond[[pond]] <- pairwise_results
}

# flatten the list into a long dataframe
fst_long <- do.call(rbind, lapply(names(fst_within_pond), function(pond) {
  mat <- fst_within_pond[[pond]]
  df <- as.data.frame(as.table(mat)) %>%
    rename(date1 = Var1, date2 = Var2, Fst = Freq) %>%
    mutate(pond = pond)
  return(df)
}))

# remove duplicate lower triangle and diagonal
fst_long <- fst_long %>%
  mutate(
    date1 = as.character(date1),
    date2 = as.character(date2)
  ) %>%
  filter(!is.na(Fst), date1 < date2)

# write to CSV

write.csv(fst_long,"/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/fst_within_pond.csv")

fst_long <- fst_long %>%
  filter(!date1 %in% c("2023", "2024"),
         !date2 %in% c("2023", "2024"))

bypondplot <- ggplot(fst_long, aes(x = date1, y = date2, fill = Fst)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.3f", Fst)), size = 3) +
  scale_fill_gradient(low = "white", high = "red", na.value = "grey90") +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    panel.grid = element_blank()
  ) +
  labs(title = "Pairwise Fst Within Ponds by Date", fill = "Fst") +
  facet_wrap(~ pond, scales = "free")

ggsave("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/popheatmap_bypondplot.png", plot = bypondplot, width = 12, height = 8, dpi = 300)















# list to hold results
fst_indiv <- list()

# loop over ponds
for (pond in unique(metadata_sub$accuratelocation)) {
  meta_p <- metadata_sub %>% filter(accuratelocation == pond)
  
  # loop over dates for this pond
  for (d in unique(meta_p$date)) {
    meta_pd <- meta_p %>% filter(date == d)
    indivs <- meta_pd$Well
    
    if (length(indivs) < 2) next  # need at least 2 individuals
    
    # initialize pairwise matrix
    pairwise_results <- matrix(NA,
                               nrow = length(indivs),
                               ncol = length(indivs),
                               dimnames = list(indivs, indivs))
    
    # loop over individual pairs
    for (i in 1:(length(indivs) - 1)) {
      for (j in (i + 1):length(indivs)) {
        
        idx <- meta_pd$Well %in% c(indivs[i], indivs[j])
        sub_samples <- meta_pd$Well[idx]
        sub_pop <- factor(meta_pd$Well[idx])  # population = individual ID
        
        if (length(unique(sub_pop)) == 2) {
          fst_tmp <- snpgdsFst(
            genofile,
            sample.id = sub_samples,
            population = sub_pop,
            snp.id = valid_variants,
            autosome.only = FALSE,
            method = "W&C84",
            maf = 0.05,
            missing.rate = 0.5,
            verbose = FALSE
          )
          
          pairwise_results[i, j] <- fst_tmp$Fst
          pairwise_results[j, i] <- fst_tmp$Fst
        }
      }
    }
    
    fst_indiv[[paste(pond, d, sep = "_")]] <- pairwise_results
  }
}





fst_long_indiv <- do.call(rbind, lapply(names(fst_indiv), function(group) {
  mat <- fst_indiv[[group]]
  if (is.null(mat)) return(NULL)
  
  parts <- strsplit(group, "_")[[1]]
  pond <- parts[1]
  date <- paste(parts[-1], collapse = "_")
  
  df <- as.data.frame(as.table(mat)) %>%
    rename(indiv1 = Var1, indiv2 = Var2, Fst = Freq) %>%
    mutate(pond = pond, date = date)
  
  return(df)
}))

fst_long_indiv <- fst_long_indiv %>%
  filter(!is.na(Fst), indiv1 != indiv2)

write.csv(fst_long_indiv,
          "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/fst_by_individual.csv",
          row.names = FALSE)


fst_long_indiv <- fst_long_indiv %>%
  mutate(indiv1 = factor(indiv1, levels = unique(indiv1)),
         indiv2 = factor(indiv2, levels = unique(indiv2)))




hist_plot <- ggplot(fst_long_indiv, aes(x = Fst)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white", alpha = 0.8) +
  facet_wrap(~ pond, scales = "free_y") +
  theme_bw() +
  labs(title = "Distribution of Pairwise Fst by Pond",
       x = "Fst", y = "Count") +
  theme(strip.text = element_text(size = 10, face = "bold"))

ggsave("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/fst_hist_by_pond_genomic.png",
       plot = hist_plot, width = 10, height = 6, dpi = 300)




fst_long_indiv <- fst_long_indiv %>%
  mutate(date = mdy(date))

fst_long_indiv <- fst_long_indiv %>%
  mutate(date = factor(date, levels = sort(unique(date))))




# Make sure date is a proper Date object
fst_long_indiv <- fst_long_indiv %>%
  mutate(date = as.character(date)) %>%  # convert to character first
  mutate(date = ifelse(nchar(date) == 4, paste0("1/1/", date), date)) %>%
  mutate(date = mdy(date)) %>%           # convert to Date
  arrange(date, pond) %>%
  mutate(
    pond = factor(pond, levels = unique(pond)),
    date = factor(date, levels = sort(unique(date)))  # chronological
  )

# Create a new column for pond + date to facet only existing combinations
fst_long_indiv <- fst_long_indiv %>%
  mutate(pond_date = paste(pond, date, sep = "_"))




fst_long_indiv <- fst_long_indiv %>%
  # Ensure date is character
  mutate(date = as.character(date)) %>%
  # Handle 4-digit years as Jan 1
  mutate(date = ifelse(grepl("^\\d{4}$", date), paste0("1/1/", date), date)) %>%
  # Parse dates safely with mdy()
  mutate(date_parsed = suppressWarnings(mdy(date))) %>%
  # If mdy fails, try ymd()
  mutate(date_parsed = ifelse(is.na(date_parsed), ymd(date), date_parsed)) %>%
  # Convert back to Date class
  mutate(date_parsed = as.Date(date_parsed, origin = "1970-01-01")) %>%
  # Create pond_type grouping
  mutate(pond_type = case_when(
    grepl("^Gilmer", pond) ~ "Gilmer",
    pond %in% c("P58") ~ "P58",
    pond %in% c("P62") ~ "P62",
    pond %in% c("P63") ~ "P63",
    pond %in% c("P66") ~ "P66",
    TRUE ~ pond
  )) %>%
  # Convert pond_type and date to factors for plotting
  mutate(
    pond_type = factor(pond_type, levels = c("Gilmer", "P58", "P62", "P63", "P66")),
    date = factor(format(date_parsed, "%m/%d/%Y"),
                  levels = unique(format(sort(date_parsed), "%m/%d/%Y")))
  )



library(ggplot2)
library(dplyr)
library(patchwork)
library(foreach)

ponds <- unique(fst_long_indiv$pond)

fst_long_indiv <- subset(fst_long_indiv, date != "01/01/2024")

pond_plots <- foreach(p = ponds, .combine = list, .multicombine = TRUE) %do% {
  
  df <- fst_long_indiv %>% filter(pond == p)
  
  ggplot(df, aes(x = indiv1, y = indiv2, fill = Fst)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                         midpoint = 0, limits = c(-1, 1), name = "Fst") +
    facet_wrap(~ date, scales = "free", nrow = 1) +  # free x and y
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      panel.grid = element_blank(),
      strip.text = element_text(size = 10, face = "bold")
    ) +
    labs(title = paste("Pond:", p),
         x = "Individual", y = "Individual")
}

# Combine vertically
final_plot <- patchwork::wrap_plots(pond_plots, ncol = 1)

# Save
ggsave("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/popheatmap_eachpond_Genomic_freexy.png",
       plot = final_plot, width = 20, height = 6 * length(ponds), dpi = 300, limitsize = FALSE)
