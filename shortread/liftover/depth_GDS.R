library(SeqArray)

# Convert VCF to GDS
#vcf.fn <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/trimmed_10bp_repeatmasked_vcf/lifted_vcf/lifted_all_missingasref.vcf.gz"          # path to your VCF file
#gds.fn <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/trimmed_10bp_repeatmasked_vcf/lifted_vcf/lifted_12major_missingasref.gds"             # output GDS file
#seqVCF2GDS(vcf.fn, gds.fn, storage.option="ZIP_RA", parallel=10, verbose=T, optimize=T)
#genofile <- seqOpen(gds.fn)

# Paths
gds.fn <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/trimmed_10bp_repeatmasked_vcf/lifted_vcf/lifted_12major_missingasref.gds"
genofile <- seqOpen(gds.fn)

# Total number of variants
n_variants <- length(seqGetData(genofile, "variant.id"))

# Subset 1 million SNPs randomly (or use head for first 1M if preferred)
set.seed(123)  # for reproducibility
subset_variants <- sample(seq_len(n_variants), size = min(1e6, n_variants))

# Set variant filter
seqSetFilter(genofile, variant.id = subset_variants, verbose = TRUE)

# Get genotype data for filtered SNPs
geno <- seqGetData(genofile, "genotype")  # [SNP x sample x ploidy]

# Missingness matrix
geno_miss <- apply(is.na(geno), c(1, 2), any)

# Per-sample missing rate
sample_missing_rate <- colMeans(geno_miss)
sample_coverage     <- colSums(!geno_miss)
sample_ids          <- seqGetData(genofile, "sample.id")

sample_df <- data.frame(
  Sample      = sample_ids,
  MissingRate = sample_missing_rate,
  Coverage    = sample_coverage
)

write.csv(
  sample_df,
  "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/trimmed_10bp_repeatmasked_vcf/lifted_vcf/per_sample_missing_coverage_1M.csv",
  row.names = FALSE
)

# Per-SNP stats
snp_missing_rate <- rowMeans(geno_miss)
snp_coverage     <- rowSums(!geno_miss)
snp_ids          <- seqGetData(genofile, "variant.id")

snp_df <- data.frame(
  VariantID   = snp_ids,
  MissingRate = snp_missing_rate,
  Coverage    = snp_coverage
)

write.csv(
  snp_df,
  "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/trimmed_10bp_repeatmasked_vcf/lifted_vcf/per_snp_missing_coverage_1M.csv",
  row.names = FALSE
)

# Close connection
seqClose(genofile)