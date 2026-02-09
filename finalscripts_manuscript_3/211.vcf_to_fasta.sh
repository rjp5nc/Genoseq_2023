#ijob -A berglandlab -c2 -p standard --mem=40G
#module load gcc/11.4.0  openmpi/4.1.4 icu R/4.3.1
#R




module load gcc/11.4.0  openmpi/4.1.4 icu R/4.3.1
R



#14642bp
library(SeqArray)
library(gdsfmt)
library(dplyr)
library(readr)
library(stringr)
library(purrr)
library(tidyr)
library(ggplot2)

gds_path <- "/scratch/rjp5nc/UK2022_2024/allsites_mito/usdobtusa_mito_allsites.annot2.ALL.gds"
mitotype_csv <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/mito_types_v2.csv"
out_dir <- "/scratch/rjp5nc/UK2022_2024/allsites_mito/unique_snps_by_group"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Mitotype map ----
mitotypes <- read.csv(mitotype_csv, stringsAsFactors = FALSE)
# Normalize columns and drop the first index column if present
stopifnot(all(c("CloneA","Group") %in% names(mitotypes)))
mito <- mitotypes %>%
  transmute(sample = as.character(CloneA),
            group  = as.character(Group)) %>%
  mutate(sample = str_trim(sample), group = str_trim(group)) %>%
  filter(nzchar(sample), nzchar(group))

# ---- Open GDS ----
options(mc.cores = 1)
Sys.setenv(OMP_NUM_THREADS = "1")

g <- seqOpen(gds_path)

## 1) Recompute per-sample mean DP (orientation-safe) and build clean keep_samples
seqResetFilter(g)
DP_all <- seqGetData(g, "annotation/format/DP"); if (is.list(DP_all)) DP_all <- do.call(cbind, DP_all)
g_samples <- seqGetData(g, "sample.id")

if (ncol(DP_all) == length(g_samples)) {
  mean_dp <- colMeans(DP_all, na.rm = TRUE)
} else if (nrow(DP_all) == length(g_samples)) {
  mean_dp <- rowMeans(DP_all, na.rm = TRUE)
} else stop("Unexpected DP matrix shape.")

mask <- (mean_dp > 40) & !is.na(mean_dp)
keep_samples <- g_samples[mask]
keep_samples <- keep_samples[!is.na(keep_samples)]  # drop any NA slots
stopifnot(all(!is.na(keep_samples)))




message(sprintf("Keeping %d/%d (mean DP>20). NA mean-DP samples: %d",
                length(keep_samples), length(g_samples), sum(is.na(mean_dp))))

## 2) Apply filter and lock the *actual* sample order
seqResetFilter(g)
seqSetFilter(g, sample.id = keep_samples, action = "set")
samps <- seqGetData(g, "sample.id")   # definitive order for columns

## 3) Rebuild mito_keep strictly from current filtered order
mitotypes <- read.csv("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/mito_types_v2.csv",
                      stringsAsFactors = FALSE)
new_samples <- c("SRR5012393", "SRR5012394", "SRR5012771",
                 "SRR5012396", "SRR5012770", "SRR5012773")

# Build a data frame to append
new_rows <- data.frame(
  X = seq(max(mitotypes$X, na.rm = TRUE) + 1,
          length.out = length(new_samples)),
  CloneA = new_samples,
  Group = "XX",
  stringsAsFactors = FALSE
)

# Append to mitotypes
mitotypes <- rbind(mitotypes, new_rows)

# Check the new tail
tail(mitotypes)

mitotypes <- subset(mitotypes, Group %in% c("A", "B", "C", "D", "E", "F", "XX"))
mitotypesE <- subset(mitotypes, Group %in% c("E", "F", "XX"))

mito_keep <- mitotypes %>%
  transmute(sample = trimws(as.character(CloneA)),
            group  = trimws(as.character(Group))) %>%
  semi_join(tibble(sample = samps), by = "sample")


write.table(mito_keep$sample, "/scratch/rjp5nc/UK2022_2024/allsites_mito/alignment_us/keep_samplesdp40.txt",
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)



write.table(mito_keep, "/scratch/rjp5nc/UK2022_2024/allsites_mito/alignment_us/conversion40.txt",
            sep = "\t",
            row.names = TRUE,
            quote = FALSE)





q()

n








module load mafft
module load samtools  # optional, only if you want quick checks

consdir="/scratch/rjp5nc/UK2022_2024/consensusmito"
alndir="/scratch/rjp5nc/UK2022_2024/allsites_mito/alignment_us"
samples="${alndir}/keep_samplesdp40.txt"

mkdir -p "$alndir"

# 1) Build list of FASTA files to include (sample + suffix)
awk '{print $0".filt.consensus.mito.fa"}' "$samples" > "${alndir}/fa_keep40.txt"

# 2) Make a multi-FASTA from those files (rename headers to sample IDs)
> "${alndir}/input_for_mafft40.fasta"
while read -r f; do
  fa="${consdir}/${f}"
  [ -f "$fa" ] || { echo "Missing: $fa" >&2; continue; }
  sample="${f%.filt.consensus.mito.fa}"
  awk -v s="$sample" '
    /^>/ {print ">" s; next}
    {gsub(/[ \t\r]/,""); if(length($0)>0) print}
  ' "$fa" >> "${alndir}/input_for_mafft40.fasta"
done < "${alndir}/fa_keep40.txt"

# 3) Align with MAFFT
mafft --auto --thread -1 "${alndir}/input_for_mafft40.fasta" > "${alndir}/alignment_aligned40.fasta"


R

library(ape)

fa <- "/scratch/rjp5nc/UK2022_2024/allsites_mito/alignment_us/alignment_aligned40.fasta"
nx <- "/scratch/rjp5nc/UK2022_2024/allsites_mito/alignment_us/alignment_aligned40.nexus"

dna <- read.dna(fa, format = "fasta")
write.nexus.data(dna, file = nx, format = "dna")

cat("Wrote", nx, "\n")





map <- read.table("/scratch/rjp5nc/UK2022_2024/allsites_mito/alignment_us/conversion40.txt")

colnames(map) <- c("old", "group")

aln <- read.nexus.data("/scratch/rjp5nc/UK2022_2024/allsites_mito/alignment_us/alignment_aligned40.nexus")


map_aln <- map[match(names(aln), map$old), ]

# Create counters within each group
map_aln$new <- ave(
  map_aln$group,
  map_aln$group,
  FUN = function(x) paste0(x, "_", seq_along(x))
)

names(aln) <- map_aln$new

write.nexus.data(aln, file = "/scratch/rjp5nc/UK2022_2024/allsites_mito/alignment_us/alignment_aligned_renamed40.nexus", format = "dna")
