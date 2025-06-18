library(SeqArray)
library(ggplot2)
library(data.table)

# Convert VCF to GDS
#vcf.fn <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/trimmed_10bp_repeatmasked_vcf/lifted_vcf/lifted_all_missingasref.vcf.gz"          # path to your VCF file
#gds.fn <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/trimmed_10bp_repeatmasked_vcf/lifted_vcf/lifted_12major_missingasref.gds"             # output GDS file
#seqVCF2GDS(vcf.fn, gds.fn, storage.option="ZIP_RA", parallel=10, verbose=T, optimize=T)
#genofile <- seqOpen(gds.fn)

gds.fn <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/trimmed_10bp_repeatmasked_vcf/lifted_vcf/lifted_12major_missingasref.gds"

genofile <- seqOpen(gds.fn)

# Get missing rates
miss_rates_all <- seqMissing(genofile, per.variant=TRUE)  
miss_rates_samps <- seqMissing(genofile, per.variant=FALSE)  

summary(miss_rates_all)
summary(miss_rates_samps)


all_var_ids <- seqGetData(genofile, "variant.id")

# Sample 100,000 random variant IDs (without replacement)
set.seed(42)  # for reproducibility


seqResetFilter(genofile, verbose = TRUE)

# Total number of variants
num_variants <- length(seqGetData(genofile, "variant.id"))

# Sample 100,000 variant indices
set.seed(42)
subset_indices <- sample(seq_len(num_variants), 100000, replace = FALSE)

# Filter by variant index (variant.sel)
variant_sel <- rep(FALSE, num_variants)
variant_sel[subset_indices] <- TRUE

# Apply the filter
seqSetFilter(genofile, variant.sel = variant_sel, verbose = TRUE)

# Confirm
sum(seqGetFilter(genofile)$variant.sel)


#missing rate for variants, subsetted to 100,000 random snps
miss_rates_all_short <- seqMissing(genofile, per.variant=TRUE)  
miss_rates_samps_short <- seqMissing(genofile, per.variant=FALSE)  

miss_rates_all_short_df <- data.frame(missing_rate = miss_rates_all_short)

miss_rates_all_short_plot <- ggplot(miss_rates_all_short_df, aes(x = missing_rate)) +
  geom_histogram(bins = 100) +
  theme_bw() +
  xlab("Missing Rate") +
  ylab("Frequency")

ggsave("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/trimmed_10bp_repeatmasked_vcf/lifted_vcf/missrates_variant.png", plot = miss_rates_all_short_plot, width = 6, height = 4, dpi = 300)



#For samples

miss_rates_samps_short_df <- data.frame(missing_rate = miss_rates_samps_short)

miss_rates_samps_short_plot <- ggplot(miss_rates_samps_short_df, aes(x = missing_rate)) +
  geom_histogram(bins = 100) +
  theme_bw() +
  xlab("Missing Rate") +
  ylab("Frequency")

# View plot

ggsave("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/trimmed_10bp_repeatmasked_vcf/lifted_vcf/miss_rates_samps_short_plot.png", plot = miss_rates_samps_short_plot, width = 6, height = 4, dpi = 300)












# Now get genotypes safely
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

write.csv(sample_df, "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/trimmed_10bp_repeatmasked_vcf/lifted_vcf/per_sample_missing_coverage_subset2.csv", row.names = FALSE)

# Per-SNP stats
snp_missing_rate <- rowMeans(geno_miss)
snp_coverage <- rowSums(!geno_miss)
snp_ids <- seqGetData(genofile, "variant.id")

snp_df <- data.frame(
  VariantID = snp_ids,
  MissingRate = snp_missing_rate,
  Coverage = snp_coverage
)

write.csv(snp_df, "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/trimmed_10bp_repeatmasked_vcf/lifted_vcf/per_snp_missing_coverage_subset2.csv", row.names = FALSE)

seqClose(genofile)



