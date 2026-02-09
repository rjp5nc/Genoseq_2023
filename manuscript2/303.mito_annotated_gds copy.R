#ijob -A berglandlab -c10 -p standard --mem=40G
#module load gcc/11.4.0  openmpi/4.1.4 icu R/4.3.1
#R

#14642bp


library(SeqArray)
library(gdsfmt)
library(ggplot2)
library(dplyr)
library(tidyr)

library(SeqArray)

# Be extra conservative with threading to avoid the fork/LZMA crash paths
options(mc.cores = 1)
Sys.setenv(OMP_NUM_THREADS = "1")

in_vcf  <- "/scratch/rjp5nc/UK2022_2024/redone_mito/euobtusa/cohort_gendb/obtusa.mito.ALLSITES.rmMissGT0.3.vcf.gz"
out_gds <- "/scratch/rjp5nc/UK2022_2024/redone_mito/euobtusa/cohort_gendb/obtusa.mito.ALLSITES.rmMissGT0.3.gds"

# Import ALL INFO fields and ALL FORMAT fields (GT + PL,DP,AD)
# storage.option=ZIP_RA is the most stable choice across builds
seqVCF2GDS(
  vcf.fn         = in_vcf,
  out.fn         = out_gds,
  storage.option = "ZIP_RA",
  parallel       = FALSE,
  info.import    = "*",                      # import all INFO keys
  fmt.import     = c("GT","DP","AD","PL"),   # all FORMAT seen in your header
  verbose        = TRUE
)

# quick sanity check
g <- seqOpen(out_gds)
print(seqSummary(g))
seqClose(g)







