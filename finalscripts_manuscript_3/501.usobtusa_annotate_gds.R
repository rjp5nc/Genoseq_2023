#ijob -A berglandlab -c10 -p standard --mem=80G
#module load gcc/11.4.0  openmpi/4.1.4 icu R/4.3.1
#R

library(SeqArray)
gds_file <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/trimmed10bp_allsites_usobtusa_renamed_annotated.gds"


seqVCF2GDS(
  vcf.fn        = "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/trimmed10bp_allsites_usobtusa_renamed_annotated.vcf.gz",
  out.fn        = "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/trimmed10bp_allsites_usobtusa_renamed_annotated.gds",
  storage.option = "ZIP_RA",
  info.import   = "ANN",
  parallel      = 8,        # ← number of threads you want
  verbose       = TRUE
)