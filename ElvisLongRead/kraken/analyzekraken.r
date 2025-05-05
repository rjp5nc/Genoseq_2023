
#module load gcc/11.4.0 openmpi/4.1.4 icu R/4.3.1
#R

library(tidyverse)

classified_headers_raw

rawreport <- read.csv("/scratch/rjp5nc/krakenDB/US_obtusa/report.txt", header=FALSE)
rawreport2 <- rawreport
rawreport2[] <- lapply(rawreport2, function(x) gsub(" ", "", x))
accessions <- separate(rawreport2, V1, into = c("percent_raw", "Numberofreads", "Numberofdirectreads", "Col4", "accession", "species"), sep = "\t")
accessions2 <- accessions[c(5,6)]

#rawheaders <- read.csv("/scratch/rjp5nc/krakenDB/US_obtusa/classified_headers_raw.txt", header=FALSE)
#hap1headers <- read.csv("/scratch/rjp5nc/krakenDB/US_obtusa/classified_headers.txt", header=FALSE)


hap2headers <- read.csv("/scratch/rjp5nc/krakenDB/US_obtusa/classified_headers.txt", header=FALSE)
hap2headers2 <- separate(hap2headers, V1, into = c("contig", "x"), sep = " ")
hap2headers3 <- separate(hap2headers2, V2, into = c("y", "accession"), sep = "\\|")


merged_df <- merge(hap2headers3, accessions2, by = "accession", all = FALSE)
onlydap <- subset(merged_df, species == "Daphniaobtusa"| species == "Daphniamagna"| species == "Daphniapulex"| species == "Daphniapulicaria")
contigs <- onlydap$contig
write(contigs, "/scratch/rjp5nc/krakenDB/US_obtusa/US_obtusa_dap_contigs.txt")


