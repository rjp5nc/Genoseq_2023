
# ijob -A berglandlab -c20 -p standard --mem=40G
### module load gcc/7.1.0  openmpi/3.1.4 R/4.1.1; R


#module load R/4.4.2;R

  library(data.table)
  library(dplyr)

euobtusa_to_eupulex <- fread("/scratch/rjp5nc/liftover/euobtusa_to_eupulex.bed")
usambigua_to_eupulex <- fread("/scratch/rjp5nc/liftover/usambigua_to_eupulex.bed")
usobtusa_to_eupulex <- fread("/scratch/rjp5nc/liftover/usobtusa_to_eupulex.bed")
uspulex_to_eupulex <- fread("/scratch/rjp5nc/liftover/uspulex_to_eupulex.bed")
eu_pulex <- fread("/scratch/rjp5nc/Reference_genomes/eu_pulex_totalHiCwithallbestgapclosed.clean.allbases.bed")


merged_df <- eu_pulex %>%
  left_join(euobtusa_to_eupulex,   by = c("V1", "V2", "V3")) %>%
  left_join(usambigua_to_eupulex,  by = c("V1", "V2", "V3")) %>%
  left_join(usobtusa_to_eupulex,   by = c("V1", "V2", "V3")) %>%
  left_join(uspulex_to_eupulex,    by = c("V1", "V2", "V3"))

write.table(
  merged_df,   # select BED columns
  file = "/scratch/rjp5nc/liftover/dictionaryoutput.bed",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)