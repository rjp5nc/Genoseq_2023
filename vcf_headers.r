
library(vcfR)

vcf_file <- "/project/berglandlab/Robert/UKSequencing2022_2024/usftp21.novogene.com/01.RawData/Bams/vcf/2022seq.concat.Removereps.vcf.gz"
vcf <- read.vcfR(vcf_file)

sample_names <- colnames(vcf@gt)[-c(1:9)]

write.csv(sample_names,"/project/berglandlab/Robert/samps.csv")

