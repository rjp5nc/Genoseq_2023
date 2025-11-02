
# ijob -A berglandlab -c20 -p standard --mem=40G
### module load gcc/7.1.0  openmpi/3.1.4 R/4.1.1; R


#module load R/4.4.2;R

#sed 's/\\t/\t/g' usambigua_to_eupulex.bed > usambigua_to_eupulex.bed


  library(data.table)
  library(dplyr)

euobtusa_to_eupulex <- fread("/scratch/rjp5nc/liftover/euobtusa_to_eupulex.bed")
usambigua_to_eupulex <- fread("/scratch/rjp5nc/liftover/usambigua_to_eupulex.bed", sep = "\t", header = FALSE, fill = TRUE)
usobtusa_to_eupulex <- fread("/scratch/rjp5nc/liftover/usobtusa_to_eupulex.bed")
uspulex_to_eupulex <- fread("/scratch/rjp5nc/liftover/uspulex_to_eupulex.bed")
eu_pulex <- fread("/scratch/rjp5nc/Reference_genomes/eu_pulex_totalHiCwithallbestgapclosed.clean.allbases.bed")





a1 <- head(euobtusa_to_eupulex, 1500)
a2 <- head(usambigua_to_eupulex, 1500)
a3 <- head(usobtusa_to_eupulex)
a4 <- head(uspulex_to_eupulex)

usambigua_to_eupulex <- subset(usambigua_to_eupulex, V4== "A"| V4== "T"| V4== "C"| V4 =="G")

merged_df <- eu_pulex %>%
  left_join(a1,   by = c("V1", "V2", "V3", "V4")) %>%
  left_join(a2,  by = c("V1", "V2", "V3", "V4")) %>%
  left_join(a3,   by = c("V1", "V2", "V3", "V4")) %>%
  left_join(a4,    by = c("V1", "V2", "V3", "V4"))

head(merged_df)


merged_df2 <- left_join(eu_pulex, euobtusa_to_eupulex, by = c("V1", "V2", "V3", "V4")) 

duplicated_rows <- ustoeu %>%
  dplyr::filter(duplicated(ustoeu[, c("V1", "V2", "V3", "V4")]) | 
                duplicated(ustoeu[, c("V1", "V2", "V3", "V4")], fromLast = TRUE))


merged_df3 <- left_join(merged_df2, usobtusa_to_eupulex, by = c("V1", "V2", "V3", "V4")) 



merged_df4 <- left_join(merged_df3, uspulex_to_eupulex, by = c("V1", "V2", "V3", "V4")) 

sum(is.na(merged_df2$V5))



ustoeu <- left_join(eu_pulex, uspulex_to_eupulex, by = c("V1", "V2", "V3", "V4")) 


  left_join(usambigua_to_eupulex,  by = c("V1", "V2", "V3", "V4")) %>%

write.table(
  merged_df2,   # select BED columns
  file = "/scratch/rjp5nc/liftover/euobtusa_to_eupulexdict.bed",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)