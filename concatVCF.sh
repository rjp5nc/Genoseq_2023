#!/usr/bin/env bash
#
#SBATCH -J concatVCF # A single job name for the array
#SBATCH --ntasks-per-node=10 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-2:00  ### 1 hours
#SBATCH --mem 24G
#SBATCH -o /scratch/rjp5nc/download_%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/download_%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

### sbatch ~/Dappupool20182019/variantCalling/concatVCF.sh
### sacct -j 22867938
### cat


module load htslib bcftools ### <- different module

 bgzip \
 -i -@ 10 \
 /project/berglandlab/Robert/UKSequencing2022_2024/usftp21.novogene.com/01.RawData/vcf/2022seq.Scaffold_1863_HRSCAF_2081.vcf \

 bgzip \
 -i -@ 10 \
  /project/berglandlab/Robert/UKSequencing2022_2024/usftp21.novogene.com/01.RawData/vcf/2022seq.Scaffold_1931_HRSCAF_2197.vcf \

 bgzip \
 -i -@ 10 \
  /project/berglandlab/Robert/UKSequencing2022_2024/usftp21.novogene.com/01.RawData/vcf/2022seq.Scaffold_2158_HRSCAF_2565.vcf \

 bgzip \
 -i -@ 10 \
  /project/berglandlab/Robert/UKSequencing2022_2024/usftp21.novogene.com/01.RawData/vcf/2022seq.Scaffold_2217_HRSCAF_2652.vcf \

 bgzip \
 -i -@ 10 \
  /project/berglandlab/Robert/UKSequencing2022_2024/usftp21.novogene.com/01.RawData/vcf/2022seq.Scaffold_2373_HRSCAF_2879.vcf \

 bgzip \
 -i -@ 10 \
  /project/berglandlab/Robert/UKSequencing2022_2024/usftp21.novogene.com/01.RawData/vcf/2022seq.Scaffold_6786_HRSCAF_7541.vcf \

 bgzip \
 -i -@ 10 \
  /project/berglandlab/Robert/UKSequencing2022_2024/usftp21.novogene.com/01.RawData/vcf/2022seq.Scaffold_7757_HRSCAF_8726.vcf \

 bgzip \
 -i -@ 10 \
  /project/berglandlab/Robert/UKSequencing2022_2024/usftp21.novogene.com/01.RawData/vcf/2022seq.Scaffold_9197_HRSCAF_10753.vcf \

 bgzip \
 -i -@ 10 \
  /project/berglandlab/Robert/UKSequencing2022_2024/usftp21.novogene.com/01.RawData/vcf/2022seq.Scaffold_9198_HRSCAF_10754.vcf \

 bgzip \
 -i -@ 10 \
  /project/berglandlab/Robert/UKSequencing2022_2024/usftp21.novogene.com/01.RawData/vcf/2022seq.Scaffold_9199_HRSCAF_10755.vcf \

 bgzip \
 -i -@ 10 \
  /project/berglandlab/Robert/UKSequencing2022_2024/usftp21.novogene.com/01.RawData/vcf/2022seq.Scaffold_9200_HRSCAF_10757.vcf \

 bgzip \
 -i -@ 10 \
  /project/berglandlab/Robert/UKSequencing2022_2024/usftp21.novogene.com/01.RawData/vcf/2022seq.Scaffold_9201_HRSCAF_10758.vcf \



 bcftools index \
  /project/berglandlab/Robert/UKSequencing2022_2024/usftp21.novogene.com/01.RawData/vcf/2022seq.Scaffold_1863_HRSCAF_2081.vcf.gz \

 bcftools index \
  /project/berglandlab/Robert/UKSequencing2022_2024/usftp21.novogene.com/01.RawData/vcf/2022seq.Scaffold_1931_HRSCAF_2197.vcf.gz \

 bcftools index \
  /project/berglandlab/Robert/UKSequencing2022_2024/usftp21.novogene.com/01.RawData/vcf/2022seq.Scaffold_2158_HRSCAF_2565.vcf.gz \

 bcftools index \
  /project/berglandlab/Robert/UKSequencing2022_2024/usftp21.novogene.com/01.RawData/vcf/2022seq.Scaffold_2217_HRSCAF_2652.vcf.gz \

 bcftools index \
  /project/berglandlab/Robert/UKSequencing2022_2024/usftp21.novogene.com/01.RawData/vcf/2022seq.Scaffold_2373_HRSCAF_2879.vcf.gz \

 bcftools index \
  /project/berglandlab/Robert/UKSequencing2022_2024/usftp21.novogene.com/01.RawData/vcf/2022seq.Scaffold_6786_HRSCAF_7541.vcf.gz \

 bcftools index \
  /project/berglandlab/Robert/UKSequencing2022_2024/usftp21.novogene.com/01.RawData/vcf/2022seq.Scaffold_7757_HRSCAF_8726.vcf.gz \

 bcftools index \
  /project/berglandlab/Robert/UKSequencing2022_2024/usftp21.novogene.com/01.RawData/vcf/2022seq.Scaffold_9197_HRSCAF_10753.vcf.gz \

 bcftools index \
  /project/berglandlab/Robert/UKSequencing2022_2024/usftp21.novogene.com/01.RawData/vcf/2022seq.Scaffold_9198_HRSCAF_10754.vcf.gz \

 bcftools index \
  /project/berglandlab/Robert/UKSequencing2022_2024/usftp21.novogene.com/01.RawData/vcf/2022seq.Scaffold_9199_HRSCAF_10755.vcf.gz \

 bcftools index \
  /project/berglandlab/Robert/UKSequencing2022_2024/usftp21.novogene.com/01.RawData/vcf/2022seq.Scaffold_9200_HRSCAF_10757.vcf.gz \

 bcftools index \
  /project/berglandlab/Robert/UKSequencing2022_2024/usftp21.novogene.com/01.RawData/vcf/2022seq.Scaffold_9201_HRSCAF_10758.vcf.gz \


bcftools concat \
--threads 10 \
-o  /project/berglandlab/Robert/UKSequencing2022_2024/usftp21.novogene.com/01.RawData/vcf/2022seq.concat.vcf.gz \
-O z \
 /project/berglandlab/Robert/UKSequencing2022_2024/usftp21.novogene.com/01.RawData/vcf/2022seq.Scaffold_1863_HRSCAF_2081.vcf.gz \
 /project/berglandlab/Robert/UKSequencing2022_2024/usftp21.novogene.com/01.RawData/vcf/2022seq.Scaffold_1931_HRSCAF_2197.vcf.gz \
 /project/berglandlab/Robert/UKSequencing2022_2024/usftp21.novogene.com/01.RawData/vcf/2022seq.Scaffold_2217_HRSCAF_2652.vcf.gz \
 /project/berglandlab/Robert/UKSequencing2022_2024/usftp21.novogene.com/01.RawData/vcf/2022seq.Scaffold_2158_HRSCAF_2565.vcf.gz \
 /project/berglandlab/Robert/UKSequencing2022_2024/usftp21.novogene.com/01.RawData/vcf/2022seq.Scaffold_2373_HRSCAF_2879.vcf.gz \
 /project/berglandlab/Robert/UKSequencing2022_2024/usftp21.novogene.com/01.RawData/vcf/2022seq.Scaffold_6786_HRSCAF_7541.vcf.gz \
 /project/berglandlab/Robert/UKSequencing2022_2024/usftp21.novogene.com/01.RawData/vcf/2022seq.Scaffold_7757_HRSCAF_8726.vcf.gz \
 /project/berglandlab/Robert/UKSequencing2022_2024/usftp21.novogene.com/01.RawData/vcf/2022seq.Scaffold_9197_HRSCAF_10753.vcf.gz \
 /project/berglandlab/Robert/UKSequencing2022_2024/usftp21.novogene.com/01.RawData/vcf/2022seq.Scaffold_9198_HRSCAF_10754.vcf.gz \
 /project/berglandlab/Robert/UKSequencing2022_2024/usftp21.novogene.com/01.RawData/vcf/2022seq.Scaffold_9199_HRSCAF_10755.vcf.gz \
 /project/berglandlab/Robert/UKSequencing2022_2024/usftp21.novogene.com/01.RawData/vcf/2022seq.Scaffold_9200_HRSCAF_10757.vcf.gz \
 /project/berglandlab/Robert/UKSequencing2022_2024/usftp21.novogene.com/01.RawData/vcf/2022seq.Scaffold_9201_HRSCAF_10758.vcf.gz
