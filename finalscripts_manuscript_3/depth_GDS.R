#module load gcc openmpi R/4.3.1


library(SeqArray)
library(ggplot2)
library(data.table)
library(tidyverse)
library(doParallel)
library(foreach)

# Convert VCF to GDS
vcf.fn <- "/scratch/rjp5nc/UK2022_2024/mito_vcf/merged_vcfs/usdobtusa_mito_genotyped.vcf.gz"          # path to your VCF file
gds.fn <- "/scratch/rjp5nc/UK2022_2024/mito_vcf/merged_vcfs/usdobtusa_mito_genotyped.gds"             # output GDS file
seqVCF2GDS(vcf.fn, gds.fn, storage.option="ZIP_RA", parallel=10, verbose=T, optimize=T)
#genofile <- seqOpen(gds.fn)

gds.fn <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/trimmed_10bp_repeatmasked_vcf/lifted_vcf/lifted_12major_missingasref.gds"

genofile <- seqOpen(gds.fn)

# Get missing rates
miss_rates_all <- seqMissing(genofile, per.variant=TRUE)  
miss_rates_samps <- seqMissing(genofile, per.variant=FALSE)  

summary(miss_rates_all)
summary(miss_rates_samps)

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

ggsave("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/trimmed_10bp_repeatmasked_vcf/lifted_vcf/missrates_variant2.png", plot = miss_rates_all_short_plot, width = 6, height = 4, dpi = 300)



#For samples

miss_rates_samps_short_df <- data.frame(missing_rate = miss_rates_samps_short)

miss_rates_samps_short_plot <- ggplot(miss_rates_samps_short_df, aes(x = missing_rate)) +
  geom_histogram(bins = 100) +
  theme_bw() +
  xlab("Missing Rate") +
  ylab("Frequency")

# View plot

ggsave("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/trimmed_10bp_repeatmasked_vcf/lifted_vcf/miss_rates_samps_short_plot2.png", plot = miss_rates_samps_short_plot, width = 6, height = 4, dpi = 300)

seqResetFilter(genofile, verbose = TRUE)


samps <- fread("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/All_clones_metadata.csv")

### get sample names
  sampleId <- seqGetData(genofile, "sample.id")

### iterate through each sample to get missing rate & coverage
  sampleStats <- foreach(sample.i=samps$Sample_ID, .combine="rbind")%dopar%{
    ## sample.i = sampleId[1]
    ### subsample to sample.i
      seqSetFilter(genofile, sample.id=sample.i)

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

  
write.csv(sampleStats, "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/trimmed_10bp_repeatmasked_vcf/lifted_vcf/sampleStats_wholedata.csv", row.names = FALSE)



#partition the samples based on species
#expect the missing rates for other species to be higher
#Realized covered , DP is depth, what fraction are 0's?
#then filter, then catagorize, build a table that has for all 8.9 mil snp's, for each of 5 species, average coverage, average depth, alt freq, classification of polymorphic, fixed sites, only looking at biallelic snp's, SeqArray there is a command to pull out sites with only 2 alleles. Build an array job for this, assign window to make it go quicker
#Package, SNPRelate works on GDS files, can do a lot of popgen analysis

#Identify clonal individuals first
#In SNPRelate, can identify identity by state. Use this to generate a heatmap
#Filter out individuals with missing rate per sample by ~5 to 10%
#Filter out variants by ~10% missing

#poppr, run within species


samps <- fread("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/All_clones_metadata.csv")
sampleId <- seqGetData(genofile, "sample.id")

# Output file path
outfile <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/trimmed_10bp_repeatmasked_vcf/lifted_vcf/sampleStats_wholedata.csv"

RobertUK_A10

# Write header
fwrite(data.table(sampleId=character(), missingRate=double(), meanDepth=double(),
                  medianDepth=double(), upper_quantile=double(), lower_quantile=double()),
       file = outfile, row.names = FALSE)

# Register parallel backend
registerDoParallel(cores = 10)  # Adjust number of cores

foreach(sample.i = sampleId, .packages = c("SeqArray", "data.table")) %dopar% {
  seqSetFilter(genofile, sample.id = sample.i)
  
  missingRate <- mean(is.na(seqGetData(genofile, "$dosage")))
  depth <- seqGetData(genofile, "annotation/format/DP")
  
  result <- data.table(
    sampleId = sample.i,
    missingRate = missingRate,
    meanDepth = mean(depth, na.rm = TRUE),
    medianDepth = median(depth, na.rm = TRUE),
    upper_quantile = quantile(depth, 0.975, na.rm = TRUE),
    lower_quantile = quantile(depth, 0.025, na.rm = TRUE)
  )
  
  # Use a lock file to prevent simultaneous writes
  lockfile <- paste0(outfile, ".lock")
  while (file.exists(lockfile)) Sys.sleep(runif(1, 0.01, 0.1))  # Wait if another process is writing
  file.create(lockfile)
  fwrite(result, file = outfile, append = TRUE, col.names = FALSE)
  file.remove(lockfile)
}

seqClose(genofile)


samps


samps <- fread("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/All_clones_metadata.csv")
samps <- fread("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/accession_plus_species.csv")
sampleId <- seqGetData(genofile, "sample.id")

missing_ids <- samps$Sample_ID[!samps$Sample_ID %in% sampleId]
missing_ids

#sampleId <- gsub("-", "_", sampleId)

missing_ids2 <- sampleId[!sampleId %in% samps$Sample_ID]
missing_ids2


missing_ids <- as.data.table(missing_ids)
missing_ids2 <- as.data.table(missing_ids2)

missingids3 <- rbind(missing_ids,missing_ids2, use.names=FALSE)

write.csv(missingids3, "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/missing.csv")
missing_ids2

write.csv(missing_ids2, "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/missing_to_put_in_csv2.csv")





unique_groups <- unique(samps$Sample_ID)

subsetsamps <- samps[,c(1,2,3,7)]

df_unique <- subsetsamps %>%
  distinct(Sample_ID, .keep_all = TRUE)


  write.csv(df_unique, "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/Sample_ID_Species.csv")


samps
samps2 <- fread("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/Sample_ID_Species.csv")


merged_data6 <- left_join(samps2, samps, by = c("Sample_ID" = "Sample_ID"))


write.csv(merged_data6, "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/Sample_ID_Species_merged.csv")

