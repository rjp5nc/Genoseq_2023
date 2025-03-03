
library(vcfR)
library(adegenet)
library(ape)



vcf <- read.vcfR("/project/berglandlab/Robert/UKSequencing2022_2024/usftp21.novogene.com/01.RawData/Bams/vcf/2022seq.concat.Removereps.vcf.gz")

# Step 2: Extract Genotype Data
genotypes <- extract.gt(vcf)

# Step 3: Convert Genotype Data to Genind Object
genind_obj <- vcfR2genind(vcf)

# Step 4: Compute Distance Matrix
dist_matrix <- dist.genpop(genind_obj)

# Step 5: Build Phylogenetic Tree (Neighbor-Joining)
tree <- nj(dist_matrix)
save(tree, file = "/project/berglandlab/Robert/UKSequencing2022_2024/phylogenetic_tree.RData")

