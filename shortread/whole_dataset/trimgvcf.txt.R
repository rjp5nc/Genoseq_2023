#module load gcc/7.1.0 openmpi/3.1.4 R/4.1.1; R

library(dplyr)

all_gvcfs <- read.csv("/scratch/rjp5nc/list_of_gvcf.csv", header=FALSE)
Allsra <- read.csv("/scratch/rjp5nc/All_clones_metadata.csv", header=TRUE)

all_gvcfs$V1 <- gsub(".g.vcf.gz", "", all_gvcfs$V1)
all_gvcfs$V1 <- gsub(".tbi", "", all_gvcfs$V1)
all_gvcfs$V1 <- gsub(".g.vcf.idx", "", all_gvcfs$V1)

all_gvcfs2 <- unique(all_gvcfs$V1)

SRA_US_pulex <- subset(Allsra, Species == "Daphnia pulex" & Continent=="NorthAmerica")

samps_already_exist_pulex <- intersect(SRA_US_pulex$Accession, all_gvcfs2)

USpulex <- filter(Allsra, Accession %in% samps_already_exist_pulex)





Euro_gvcfs <- read.csv("/scratch/rjp5nc/Euro_gvcfs.csv", header=FALSE)

Euro_gvcfs$V1 <- gsub(".g.vcf.gz", "", Euro_gvcfs$V1)
Euro_gvcfs$V1 <- gsub(".tbi", "", Euro_gvcfs$V1)
Euro_gvcfs$V1 <- gsub(".g.vcf.idx", "", Euro_gvcfs$V1)

Euro_gvcfs2 <- unique(Euro_gvcfs$V1)

SRA_Eu_pulex <- subset(Allsra, Species == "Daphnia pulex" & Continent == "Europe")

samps_already_exist_Eupulex <- intersect(SRA_Eu_pulex$Sample_ID, Euro_gvcfs2)


Eupulex <- filter(Allsra, Sample_ID %in% samps_already_exist_Eupulex)


existing_gvcf <- rbind(Eupulex,USpulex)

need <- Allsra %>% filter(!Sample_ID %in% existing_gvcf$Sample_ID)

forref2 <- subset(need, Accession != "")

downloadsra2 <- downloadsra$Accession

write.csv(downloadsra2,"/scratch/rjp5nc/downloadsra.csv")
write.csv(forref2,"/scratch/rjp5nc/forref2.csv")
