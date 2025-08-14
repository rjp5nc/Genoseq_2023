#module load gcc openmpi R/4.3.1
#R

library(SeqArray)
library(SeqVarTools)
library(ggplot2)
library(doMC)
library(data.table)

gds.fn <- "/scratch/rjp5nc/UK2022_2024/buscoanalysis/buscolifted_singlecopy.gds.gz"

samples <- read.csv("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/Sample_ID_Species_merged_20250627.csv")

genofile <- seqOpen(gds.fn)

### Depth filter step — keep only variants where all samples have mean depth > 1
variant_depth <- seqApply(genofile, "annotation/format/DP", 
                          margin = "by.variant", 
                          FUN = function(dp) mean(dp, na.rm = TRUE), 
                          as.is = "double")

variant_depth_filter <- variant_depth > 1
cat("Variants passing depth > 1 filter:", sum(variant_depth_filter), "\n")

# Apply depth filter to GDS
seqSetFilter(genofile, variant.sel = variant_depth_filter, verbose = TRUE)

### Get missing rates per sample after depth filter
miss_rates_samps <- seqMissing(genofile, per.variant = FALSE)
summary(miss_rates_samps)

### Get total number of variants after filtering
num_variants <- sum(variant_depth_filter)
cat("Number of variants after depth filter:", num_variants, "\n")

# Randomly subset to 100,000 variants (or all if fewer)
set.seed(22)
subset_indices <- sample(which(variant_depth_filter), 
                         min(100000, num_variants), 
                         replace = FALSE)

variant_sel <- rep(FALSE, length(variant_depth_filter))
variant_sel[subset_indices] <- TRUE

seqSetFilter(genofile, variant.sel = variant_sel, verbose = TRUE)

### Missing rates for variants (depth > 1 and subsetted)
miss_rates_all_short <- seqMissing(genofile, per.variant = TRUE)
miss_rates_all_short_df <- data.frame(missing_rate = miss_rates_all_short)

miss_rates_all_short_plot <- ggplot(miss_rates_all_short_df, aes(x = missing_rate)) +
  geom_histogram(bins = 100) +
  theme_bw() + 
  xlab("Missing Rate") +
  ylab("Frequency")

ggsave("/scratch/rjp5nc/UK2022_2024/buscoanalysis/busco_missrates_variant_depth>1.png", 
       plot = miss_rates_all_short_plot, width = 6, height = 4, dpi = 300)

### Missing rates per sample
miss_rates_samps_short <- seqMissing(genofile, per.variant = FALSE)
miss_rates_samps_short_df <- data.frame(missing_rate = miss_rates_samps_short)

miss_rates_samps_short_plot <- ggplot(miss_rates_samps_short_df, aes(x = missing_rate)) +
  geom_histogram(bins = 100) +
  theme_bw() +
  xlab("Missing Rate") +
  ylab("Frequency")

ggsave("/scratch/rjp5nc/UK2022_2024/buscoanalysis/missrates_samps_depth>1.png", 
       plot = miss_rates_samps_short_plot, width = 6, height = 4, dpi = 300)

### Per-sample depth and missing stats
seqResetFilter(genofile, verbose = TRUE)
seqSetFilter(genofile, variant.sel = variant_depth_filter, verbose = TRUE)

registerDoMC(10)

sampleId <- seqGetData(genofile, "sample.id")

sampleStats <- foreach(sample.i = sampleId, .combine = "rbind") %dopar% {
  message(sample.i)
  seqSetFilter(genofile, variant.sel = variant_depth_filter, sample.id = sample.i, verbose = FALSE)

  missingRate <- mean(is.na(seqGetData(genofile, "$dosage")))

  depth <- seqGetData(genofile, "annotation/format/DP")

  data.table(
    sampleId = sample.i,
    missingRate = missingRate,
    meanDepth = mean(depth, na.rm = TRUE),
    medianDepth = median(depth, na.rm = TRUE),
    upper_quantile = quantile(depth, .975, na.rm = TRUE),
    lower_quantile = quantile(depth, .025, na.rm = TRUE)
  )
}

write.csv(sampleStats, "/scratch/rjp5nc/UK2022_2024/buscoanalysis/sampleStats_busco_depth>1.csv", row.names = FALSE)