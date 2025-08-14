#module load gcc/7.1.0 openmpi/3.1.4 R/4.1.1; R

library(SeqArray)


#seqVCF2GDS("/project/berglandlab/Robert/UKSequencing2022_2024/usftp21.novogene.com/01.RawData/Bams/vcf/2022seq.concat.Removereps.vcf.gz", "/project/berglandlab/Robert/UKSequencing2022_2024/usftp21.novogene.com/01.RawData/Bams/vcf/2022seq.gds", storage.option="ZIP_RA")

args = commandArgs(trailingOnly=TRUE)
vcf.fn=args[[1]]
gds.fn=gsub(".vcf", ".gds", vcf.fn)

output_dir = "/scratch/rjp5nc/UK2022_2024/buscoanalysis/"
gds.fn = paste(output_dir, basename(gds.fn), sep="")

#vcf.fn=paste(vcf.fn, ".gz", sep="")
#vcf.fn="/project/berglandlab/Robert/UKSequencing2022_2024/usftp21.novogene.com/01.RawData/Bams/vcf/2022seq.concat.Removereps.vcf.gz"
seqParallelSetup(cluster=10, verbose=TRUE)

seqVCF2GDS(vcf.fn, gds.fn, storage.option="ZIP_RA", parallel=10, verbose=T, optimize=T)
