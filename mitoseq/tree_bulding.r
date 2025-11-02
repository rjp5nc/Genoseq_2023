
#module load gcc/11.4.0  openmpi/4.1.4 icu R/4.3.1


library(vcfR)
library(adegenet)
library(ape)
library(dplyr)
library(proxy)

vcf <- read.vcfR("/scratch/rjp5nc/mitobams/mito.Scaffold_7757_HRSCAF_8726.vcf")

# Step 2: Extract Genotype Data
genotypes <- extract.gt(vcf)

# Step 3: Convert Genotype Data to Genind Object
genind_obj <- vcfR2genind(vcf)


metadata <- read.csv("/project/berglandlab/Robert/UKSequencing2022_2024/2022_2024seqmetadata20250131.csv", header = TRUE)
metadata_with_clone <- read.csv("/project/berglandlab/Robert/UKSequencing2022_2024/metadata_with_clone.csv", header = TRUE)

head(metadata)

samples_to_keep <- metadata_with_clone %>% filter(accuratelocation == "P63"|accuratelocation == "P62"|accuratelocation == "PBO66"|accuratelocation == "P58") %>% pull(Sample)
  
#samples_to_keep <- metadata_with_clone %>% pull(Sample)
#
#|accuratelocation == "P63"|accuratelocation == "P62"|accuratelocation == "PBO66"|accuratelocation == "P58"|clone=="simo"
# & Plate != "simo"

## Step 6: Subset Genind Object to Keep Only Selected Samples
genind_obj2 <- genind_obj[indNames(genind_obj) %in% samples_to_keep, ]

gen_matrix <- tab(genind_obj2)

# Step 4: Compute Distance Matrix

dist_matrix <- dist(gen_matrix)


# Step 5: Build Phylogenetic Tree (Neighbor-Joining)
treemito <- njs(dist_matrix)
save(treemito, file = "/project/berglandlab/Robert/UKSequencing2022_2024/phylogenetic_treemitoobtusa.RData")


#sim_matrix <- 1 / (1 + (dist_matrix))  # Convert to similarity (scaled 0-1)

#genind_obj3<- tab(genind_obj2)
#sim_matrix <- as.matrix(simil(genind_obj, method = "Jaccard"))

heatmap(sim_matrix)  # Quick visualization




genind_data <- tab(genind_obj2)  # Convert genind object to genotype table
sim_matrix <- simil(genind_data, method = "Jaccard")

# Optionally, convert to a numeric matrix for further processing
sim_matrix <- as.matrix(sim_matrix)

write.csv(sim_matrix, "/project/berglandlab/Robert/UKSequencing2022_2024/similarity_matrixobtusa.csv")  # Save as CSV
