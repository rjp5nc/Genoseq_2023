library(SeqArray)

# Convert VCF to GDS
vcf.fn <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/trimmed_10bp_repeatmasked_vcf/lifted_vcf/lifted_all_missingasref.vcf.gz"          # path to your VCF file
gds.fn <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/trimmed_10bp_repeatmasked_vcf/lifted_vcf/lifted_12major_missingasref.gds"             # output GDS file

seqVCF2GDS(vcf.fn, gds.fn, storage.option="ZIP_RA", parallel=10, verbose=T, optimize=T)



genofile <- seqOpen(gds.fn)


geno <- seqGetData(genofile, "genotype")  # 3D array: SNP x sample x ploidy

# Missingness: TRUE if any ploidy is NA
geno_miss <- apply(is.na(geno), c(1,2), any)

# Per-sample missing rate
sample_missing_rate <- colMeans(geno_miss)

# Coverage: number of non-missing genotypes per sample
sample_coverage <- colSums(!geno_miss)

# Sample names
sample_ids <- seqGetData(genofile, "sample.id")

# Combine into data frame
sample_df <- data.frame(
  Sample = sample_ids,
  MissingRate = sample_missing_rate,
  Coverage = sample_coverage
)

# Write to CSV
write.csv(sample_df, "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/trimmed_10bp_repeatmasked_vcf/lifted_vcf/per_sample_missing_coverage.csv", row.names = FALSE)

# Per-SNP missing rate
snp_missing_rate <- rowMeans(geno_miss)

# Coverage: number of non-missing genotypes per SNP
snp_coverage <- rowSums(!geno_miss)

# SNP IDs (can also use pos/chrom)
snp_ids <- seqGetData(genofile, "variant.id")

# Combine into data frame
snp_df <- data.frame(
  VariantID = snp_ids,
  MissingRate = snp_missing_rate,
  Coverage = snp_coverage
)

# Write to CSV
write.csv(snp_df, "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/trimmed_10bp_repeatmasked_vcf/lifted_vcf/per_snp_missing_coverage.csv", row.names = FALSE)

seqClose(genofile)
