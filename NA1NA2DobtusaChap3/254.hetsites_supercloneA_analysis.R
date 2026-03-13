#ijob -A berglandlab -c2 -p standard --mem=80G
#module load gcc/11.4.0  openmpi/4.1.4 icu R/4.3.1
#R
library(SeqArray)

library(data.table)
library(ggplot2)
library(ggtern)

### paths
gds_path <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/trimmed10bp_filtered_Gilmer_12.subset_regions.gds"
meta_path <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/genomic_types.csv"
outdir    <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/"

genofile <- seqOpen(gds_path)

## --- A-only sample list (your CSV has a leading header) ---
geno_types <- fread(meta_path, skip = 1, header = FALSE,
                    col.names = c("idx","CloneA","Group"))
A_ids <- geno_types[Group == "A", CloneA]

## --- filter to A samples present in this GDS ---
sids     <- seqGetData(genofile, "sample.id")
A_in_gds <- intersect(sids, A_ids)
stopifnot(length(A_in_gds) > 0)

seqResetFilter(genofile)
seqSetFilter(genofile, sample.id = A_in_gds)




# Total number of variants
num_variants <- length(seqGetData(genofile, "variant.id"))
num_variants
# Sample 100,000 variant indices
set.seed(22)
subset_indices <- sample(seq_len(num_variants), 20000, replace = FALSE)

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


