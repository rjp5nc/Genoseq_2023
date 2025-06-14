library(SeqArray)

# Convert VCF to GDS
#vcf.fn <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/trimmed_10bp_repeatmasked_vcf/lifted_vcf/lifted_all_missingasref.vcf.gz"          # path to your VCF file
#gds.fn <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/trimmed_10bp_repeatmasked_vcf/lifted_vcf/lifted_12major_missingasref.gds"             # output GDS file
#seqVCF2GDS(vcf.fn, gds.fn, storage.option="ZIP_RA", parallel=10, verbose=T, optimize=T)
#genofile <- seqOpen(gds.fn)

gds.fn <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/trimmed_10bp_repeatmasked_vcf/lifted_vcf/lifted_12major_missingasref.gds"

# Open GDS
genofile <- seqOpen(gds.fn)

# Total number of variants
#total_variants <- seqSummary(genofile)$genotype$dim[1]

# Sample 1 million variant IDs
set.seed(42)  # for reproducibility


miss_rates_all <- seqMissing(genofile)  # returns % missing per SNP
keep_variants <- which(miss_rates_all < 1)
subset_variants <- sample(keep_variants, size = 1e6)

#subset_variants <- sample(seq_len(total_variants), size = 1e6)

# Apply SNP filter
seqSetFilter(genofile, variant.id = subset_variants, verbose = TRUE)

# Get genotype array (SNP x Sample x Ploidy)
geno <- seqGetData(genofile, "genotype")

# Compute missingness (TRUE if any ploidy is NA)
geno_miss <- apply(is.na(geno), c(1,2), any)

# Per-sample missing rate and coverage
sample_missing_rate <- colMeans(geno_miss)
sample_coverage <- colSums(!geno_miss)
sample_ids <- seqGetData(genofile, "sample.id")

sample_df <- data.frame(
  Sample = sample_ids,
  MissingRate = sample_missing_rate,
  Coverage = sample_coverage
)

write.csv(sample_df, "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/trimmed_10bp_repeatmasked_vcf/lifted_vcf/per_sample_missing_coverage_subset.csv", row.names = FALSE)

# Per-SNP missing rate and coverage
snp_missing_rate <- rowMeans(geno_miss)
snp_coverage <- rowSums(!geno_miss)
snp_ids <- seqGetData(genofile, "variant.id")

snp_df <- data.frame(
  VariantID = snp_ids,
  MissingRate = snp_missing_rate,
  Coverage = snp_coverage
)

write.csv(snp_df, "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/trimmed_10bp_repeatmasked_vcf/lifted_vcf/per_snp_missing_coverage_subset.csv", row.names = FALSE)

# Done
seqClose(genofile)
