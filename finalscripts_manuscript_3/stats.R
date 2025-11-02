#module load gcc openmpi R/4.3.1
#R

library(SeqArray)
library(ggplot2)
library(tidyverse)
library(doMC)
library(foreach)
library(data.table)


vcf.fn <- "/project/berglandlab/Robert/UKSequencing2022_2024/raw_vcfs/raw_usambigua_vcf.vcf.gz"          # path to your VCF file
gds.fn <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/trimmed_10bp_repeatmasked_vcf/raw_usambigua.gds"             # output GDS file
seqVCF2GDS(vcf.fn, gds.fn, storage.option="ZIP_RA", parallel=10, verbose=T, optimize=T)



#gds.fn <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/trimmed_10bp_repeatmasked_vcf/lifted_vcf/lifted_12major_missingasref.gds"
gds.fn <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/trimmed_10bp_repeatmasked_vcf/trimmed10bp_masked_usambigua.gds"
samples <- read.csv("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/Sample_ID_Species_merged_20250627.csv")

genofile <- seqOpen(gds.fn)

Species <- subset(samples, Species == "Daphnia ambigua" & Continent =="NorthAmerica")
Species_Id <- Species$Sample_ID

seqResetFilter(genofile, verbose = TRUE)
seqSetFilter(genofile, sample.id = Species_Id, verbose = TRUE)

Species_Id
# Get missing rates
miss_rates_samps <- seqMissing(genofile, per.variant=FALSE)  

summary(miss_rates_samps)

# Total number of variants
num_variants <- length(seqGetData(genofile, "variant.id"))
num_variants
# Sample 100,000 variant indices
set.seed(22)
subset_indices <- sample(seq_len(num_variants), 100000, replace = FALSE)

# Filter by variant index (variant.sel)
variant_sel <- rep(FALSE, num_variants)
variant_sel[subset_indices] <- TRUE

# Apply the filter
seqSetFilter(genofile, variant.sel = variant_sel, sample.id=Species_Id, verbose = TRUE)

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

ggsave("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/missrates_variant_ambigua_raw_beforeeverything.png", plot = miss_rates_all_short_plot, width = 6, height = 4, dpi = 300)



#For samples

miss_rates_samps_short_df <- data.frame(missing_rate = miss_rates_samps_short)

miss_rates_samps_short_plot <- ggplot(miss_rates_samps_short_df, aes(x = missing_rate)) +
  geom_histogram(bins = 100) +
  theme_bw() +
  xlab("Missing Rate") +
  ylab("Frequency")

# View plot

ggsave("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/missrates_variant_USambigua_raw_samps_beforeall.png", plot = miss_rates_samps_short_plot, width = 6, height = 4, dpi = 300)

seqResetFilter(genofile, verbose = TRUE)

seqSetFilter(genofile, variant.sel = variant_sel, sample.id=Species_Id, verbose = TRUE)


registerDoMC(10)


### get sample names
  sampleId <- seqGetData(genofile, "sample.id")

### iterate through each sample to get missing rate & coverage

sampleStats <- foreach(sample.i=sampleId, .combine = "rbind")%dopar%{
      # samp.i=sampleStats$sampleId[1]
      message(sample.i)
      seqSetFilter(genofile, variant.sel = variant_sel, sample.id=sample.i, verbose = TRUE)

    ### get missing data rate
      missingRate <- mean(is.na(seqGetData(genofile, "$dosage")))

    ### coverage
      depth <- seqGetData(genofile, "annotation/format/DP")

    ### output
      data.table(sampleId=sample.i,
                 missingRate=missingRate,
                 meanDepth=mean(depth, na.rm=T),
                 medianDepth=median(depth, na.rm=T),
                 upper_quantile=quantile(depth, .975, na.rm=T),
                 lower_quantile=quantile(depth, .025, na.rm=T))
  }

  
write.csv(sampleStats, "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/sampleStats_ambigua_raw_beforeeverything.csv", row.names = FALSE)

