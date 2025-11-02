#ijob -A berglandlab -c10 -p standard --mem=40G
#module load gcc/11.4.0  openmpi/4.1.4 icu R/4.3.1
#R
library(data.table)
pops <- read.table("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/pops_fixed.txt",
                   header = FALSE, stringsAsFactors = FALSE)

genomic_types <-   read.csv("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/genomic_types.csv")


genomic_ids <- genomic_types$CloneA
pop_ids <- pops$V1

# find those in pops but not in genomic_types
pops_not_in_genomic <- setdiff(pop_ids, genomic_ids)

# make a dataframe for these, with Group = "OO"
new_rows <- data.frame(
  CloneA = pops_not_in_genomic,
  Group = "OO",
  stringsAsFactors = FALSE
)

# bind to genomic_types
genomic_types_extended <- rbind(
  genomic_types[, c("CloneA", "Group")],  # keep only relevant cols
  new_rows
)

# if you want to keep the "X" column as well (row numbers), you can reset like this:
genomic_types_extended$X <- seq_len(nrow(genomic_types_extended))

# reorder columns to match original
genomic_types_extended <- genomic_types_extended[, c("X", "CloneA", "Group")]



depths <-   read.csv("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/sampleStats_US_obtusa.csv")

library(dplyr)

genomic_types_final <- genomic_types_extended %>%
  left_join(depths, by = c("CloneA" = "sampleId"))


mitotypes <-   read.csv("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/mito_types.csv")

genomic_types_final2 <- genomic_types_final %>%
  left_join(mitotypes, by = c("CloneA" = "CloneA"))

# make sure it's a data.table
setDT(genomic_types_final2)

genomic_types_final <- subset(genomic_types_final2, Group.y !="NA")

setnames(genomic_types_final, "Group.x", "Genomic_type")
setnames(genomic_types_final, "Group.y", "Mito_type")
setDT(genomic_types_final)

# split out OO's
oo_group <- genomic_types_final[Genomic_type == "OO"]

# for all other groups: order by meanDepth descending, then take top 2
top_by_group <- genomic_types_final[Genomic_type != "OO",
                                    .SD[order(-meanDepth)][1],
                                    by = Genomic_type]

top_by_group <- na.omit(top_by_group)

# combine results
final_selection <- rbind(top_by_group, oo_group)

write.csv(final_selection, "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/one_of_each_genomic_type.csv")


final_vcf_dt <- data.table(final_selection$CloneA)

# write to file
fwrite(final_vcf_dt, 
       "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/final_vcf_filter_one_of_each.txt",
       sep = "\t", quote = FALSE, row.names = FALSE)
