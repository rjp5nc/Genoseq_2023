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
library(doParallel)
registerDoParallel(8)


library(patchwork)
library(foreach)
library(lubridate)


# ---- Step 1: Open GDS file ----

metadata <- read.csv("/project/berglandlab/Robert/UKSequencing2022_2024/old_stuff/2022_2024seqmetadata20250131.csv", header = TRUE)
metadata_with_clone <- read.csv("/project/berglandlab/Robert/UKSequencing2022_2024/old_stuff/metadata_with_clone.csv", header = TRUE)
samples <- read.csv("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/Sample_ID_Species_merged_20250627.csv")
usobtusasamps <- subset(samples, Species == "Daphnia obtusa" & Continent == "NorthAmerica")
mitotypes <- read.csv("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/mito_types.csv")

metadata_with_clone <- subset(metadata_with_clone, clone !="Blank")
metadata_with_clone <- subset(metadata_with_clone, clone !="BLANK")

head(metadata)
metadata$clone <- trimws(metadata$clone)
metadata <- metadata %>% 
  filter(!tolower(clone) %in% c("blank", "blanks", "na", "missing"))


#seqClose(genofile)  
genofile.fn <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/trimmed10bp_masked_usobtusa.gds"
genofile <- seqOpen(genofile.fn)


seqResetFilter(genofile)

samples_to_keep <- usobtusasamps %>% filter(Species == "Daphnia obtusa") %>% pull(Sample_ID)
#samples_to_keep <- metadata_with_clone %>% filter(accuratelocation == "P58") %>% pull(Well)
#samples_to_keep <- metadata_with_clone %>% filter(location == "UK" | accuratelocation == "P759") %>% pull(Well)

#samples_to_keep <- metadata_with_clone$Well

unique(metadata_with_clone$accuratelocation)

seqSetFilter(genofile, sample.id = samples_to_keep)

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
write.csv(dist_matrix, "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/usobtusa_similarity_matrix2.csv")  # Save as CSV

dist_matrix2 <- 1- dist_matrix



long_mat <- melt(dist_matrix2, varnames = c("CloneA", "CloneB"), value.name = "Similarity") %>%
  mutate(CloneA = as.character(CloneA),
         CloneB = as.character(CloneB))  # ensure they are character, not factor

# 2. Filter dataset
long_filt <- long_mat %>%
  filter(Similarity >= 0.95, Similarity != 1)

g <- graph_from_data_frame(long_filt[, c("CloneA", "CloneB")], directed = FALSE)

# 4. Find connected components
comp <- components(g)

# 5. Assign group letters
group_letters <- setNames(LETTERS[comp$membership], names(comp$membership))

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

metadata_sub <- metadata %>%
  filter(Well %in% final_valid_samples)
metadata_sub <- metadata_sub[metadata_sub$accuratelocation != "PBO66", ]
metadata_sub <- metadata_sub[metadata_sub$accuratelocation != "", ]

unique(metadata_sub$accuratelocation)

pop_factor <- factor(metadata_sub$accuratelocation)

seqResetFilter(genofile)

seqSetFilter(genofile, sample.id = metadata_sub$Well, variant.id = valid_variants)

write.csv(metadata_sub, "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/metadata_sub2.csv")  # Save as CSV

metadata_sub_sum <- metadata_sub %>%
  group_by(accuratelocation, date) %>%
  summarise(count = n(), .groups = "drop")


#FST overall between ponds

pop_levels <- unique(metadata_sub$accuratelocation)
pairwise_results <- matrix(NA, nrow=length(pop_levels), ncol=length(pop_levels),
                           dimnames=list(pop_levels, pop_levels))

for(i in 1:(length(pop_levels)-1)){
  for(j in (i+1):length(pop_levels)){
    # logical index for samples in populations i and j

      message("Processing: Pop=", pop_levels)

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

write.csv(pairwise_results,"/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/pairwise_fst2.csv")

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
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-1, 1),
                       name = "Fst") +  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    panel.grid = element_blank()
  ) +
  labs(title = "Pairwise Fst Heatmap", fill = "Fst")


ggsave("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/popheatmap_genomic2.png", plot = popheatmapplot, width = 7, height = 6, dpi = 300)












#FST by each date

library(dplyr)


metadata_sub <- metadata %>%
  filter(Well %in% final_valid_samples)
metadata_sub <- metadata_sub[metadata_sub$accuratelocation != "PBO66", ]
metadata_sub <- metadata_sub[metadata_sub$accuratelocation != "", ]

#metadata_sub <- subset(metadata_sub, accuratelocation == "P62"| accuratelocation == "P66" )

unique(metadata_sub$accuratelocation)

seqResetFilter(genofile)

seqSetFilter(genofile, sample.id = metadata_sub$Well, variant.id = valid_variants)

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
  filter(!is.na(Fst))

# write to CSV
write.csv(fst_long_date,
          "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/fst_by_date_genomic.csv",
          row.names = FALSE)


bydateplot <- ggplot(fst_long_date, aes(x = pop1, y = pop2, fill = Fst)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.3f", Fst)), size = 3) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-1, 1),
                       name = "Fst") +  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    panel.grid = element_blank()
  ) +
  labs(title = "Pairwise Fst by Date", fill = "Fst") +
  facet_wrap(~ date)


ggsave("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/popheatmap_bydateplot_genomic2.png", plot = bydateplot, width = 7, height = 6, dpi = 300)






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

      message("Processing: Pond=", pond,
                " | Date=", date_levels)
      
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
  filter(!is.na(Fst))

# write to CSV

write.csv(fst_long,"/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/fst_within_pond_genomic.csv")

fst_long <- fst_long %>%
  filter(!date1 %in% c("2023", "2024"),
         !date2 %in% c("2023", "2024"))

bypondplot <- ggplot(fst_long, aes(x = date1, y = date2, fill = Fst)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.3f", Fst)), size = 3) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-1, 1),
                       name = "Fst") +  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    panel.grid = element_blank()
  ) +
  labs(title = "Pairwise Fst Within Ponds by Date", fill = "Fst") +
  facet_wrap(~ pond, scales = "free")

ggsave("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/popheatmap_bypondplot_genomic2.png", plot = bypondplot, width = 12, height = 8, dpi = 300)















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
        
        # ---- print progress ----
        message("Processing: Pond=", pond,
                " | Date=", d,
                " | ", indivs[i], " vs ", indivs[j])
        
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
  filter(!is.na(Fst))

write.csv(fst_long_indiv,
          "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/fst_by_individual_genomic2.csv",
          row.names = FALSE)
fst_long_indiv <- read.csv("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/fst_by_individual_genomic2.csv")

fst_long_indiv <- fst_long_indiv %>%
  mutate(indiv1 = factor(indiv1, levels = unique(indiv1)),
         indiv2 = factor(indiv2, levels = unique(indiv2)))

hist_plot <- ggplot(fst_long_indiv, aes(x = Fst)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white", alpha = 0.8) +
  facet_wrap(~ pond) + xlim(-1,1)+ 
  theme_bw() +
  labs(title = "Distribution of Pairwise Fst by Pond",
       x = "Fst", y = "Count") +
  theme(strip.text = element_text(size = 10, face = "bold"))

ggsave("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/fst_hist_by_pond_genomic2.png",
       plot = hist_plot, width = 10, height = 6, dpi = 300)



fst_long_indiv <- fst_long_indiv %>%
  mutate(date = as.Date(date, format = "%m/%d/%Y")) %>%  # back to Date for ordering
  arrange(date) %>%
  mutate(date = factor(format(date, "%m/%d/%Y"), 
                       levels = unique(format(date, "%m/%d/%Y"))))



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


ponds <- unique(fst_long_indiv$pond_type)

fst_long_indiv <- subset(fst_long_indiv, date != "01/01/2024")

# Run foreach in parallel
pond_plots <- foreach(p = ponds, .packages = c("ggplot2", "dplyr")) %dopar% {
  
  df <- fst_long_indiv %>% filter(pond_type == p)
  
  ggplot(df, aes(x = indiv1, y = indiv2, fill = Fst)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                         midpoint = 0, limits = c(-1, 1), name = "Fst") +
    facet_wrap(~ date, scales = "free", nrow = 1) +
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
ggsave(
  "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/popheatmap_eachpond_Genomic_freexy2.png",
  plot = final_plot, width = 20, height = 6 * length(ponds), dpi = 300, limitsize = FALSE
)




summary(subset(fst_long_indiv, pond == "P66")$Fst)















