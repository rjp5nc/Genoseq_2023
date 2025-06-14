library(SeqArray)

# Convert VCF to GDS
#vcf.fn <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/trimmed_10bp_repeatmasked_vcf/lifted_vcf/lifted_all_missingasref.vcf.gz"          # path to your VCF file
#gds.fn <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/trimmed_10bp_repeatmasked_vcf/lifted_vcf/lifted_12major_missingasref.gds"             # output GDS file
#seqVCF2GDS(vcf.fn, gds.fn, storage.option="ZIP_RA", parallel=10, verbose=T, optimize=T)
#genofile <- seqOpen(gds.fn)

gds.fn <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/trimmed_10bp_repeatmasked_vcf/lifted_vcf/lifted_12major_missingasref.gds"

genofile <- seqOpen(gds.fn)

# Get variant-level missing rates (returns proportions)
miss_rates_all <- seqMissing(genofile, per.variant=TRUE)  

# Filter out fully missing variants
keep_variants <- which(miss_rates_all < 1)

# Double-check that there are enough variants
if (length(keep_variants) < 1e6) {
  stop("Not enough non-missing SNPs to sample 1 million.")
}

# Randomly sample 1 million SNPs from non-missing
set.seed(42)
subset_variants <- sample(keep_variants, size = 1e6)

# Apply filter **before** loading genotypes
seqSetFilter(genofile, variant.id = subset_variants, verbose = TRUE)

# Get genotype array (SNP x Sample x Ploidy)
geno <- seqGetData(genofile, "genotype")

# Missingness logic (any NA across ploidy)
geno_miss <- apply(is.na(geno), c(1, 2), any)

# Per-sample stats
sample_missing_rate <- colMeans(geno_miss)
sample_coverage <- colSums(!geno_miss)
sample_ids <- seqGetData(genofile, "sample.id")

sample_df <- data.frame(
  Sample = sample_ids,
  MissingRate = sample_missing_rate,
  Coverage = sample_coverage
)

write.csv(sample_df, "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/trimmed_10bp_repeatmasked_vcf/lifted_vcf/per_sample_missing_coverage_subset.csv", row.names = FALSE)

# Per-SNP stats
snp_missing_rate <- rowMeans(geno_miss)
snp_coverage <- rowSums(!geno_miss)
snp_ids <- seqGetData(genofile, "variant.id")

snp_df <- data.frame(
  VariantID = snp_ids,
  MissingRate = snp_missing_rate,
  Coverage = snp_coverage
)

write.csv(snp_df, "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/trimmed_10bp_repeatmasked_vcf/lifted_vcf/per_snp_missing_coverage_subset.csv", row.names = FALSE)

seqClose(genofile)



