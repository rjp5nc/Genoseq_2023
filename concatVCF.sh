#!/usr/bin/env bash
#
#SBATCH -J concatVCF # A single job name for the array
#SBATCH --ntasks-per-node=10 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-4:00  ### 1 hours
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
 /scratch/rjp5nc/UK2022_2024/allshortreads/2022seq.Scaffold_1863_HRSCAF_2081.renamed.vcf \

 bgzip \
 -i -@ 10 \
  /scratch/rjp5nc/UK2022_2024/allshortreads/2022seq.Scaffold_1931_HRSCAF_2197.renamed.vcf \

 bgzip \
 -i -@ 10 \
  /scratch/rjp5nc/UK2022_2024/allshortreads/2022seq.Scaffold_2158_HRSCAF_2565.renamed.vcf \

 bgzip \
 -i -@ 10 \
  /scratch/rjp5nc/UK2022_2024/allshortreads/2022seq.Scaffold_2217_HRSCAF_2652.renamed.vcf \

 bgzip \
 -i -@ 10 \
  /scratch/rjp5nc/UK2022_2024/allshortreads/2022seq.Scaffold_2373_HRSCAF_2879.renamed.vcf \

 bgzip \
 -i -@ 10 \
  /scratch/rjp5nc/UK2022_2024/allshortreads/2022seq.Scaffold_6786_HRSCAF_7541.renamed.vcf \

 bgzip \
 -i -@ 10 \
  /scratch/rjp5nc/UK2022_2024/allshortreads/2022seq.Scaffold_7757_HRSCAF_8726.renamed.vcf \

 bgzip \
 -i -@ 10 \
  /scratch/rjp5nc/UK2022_2024/allshortreads/2022seq.Scaffold_9197_HRSCAF_10753.renamed.vcf \

 bgzip \
 -i -@ 10 \
  /scratch/rjp5nc/UK2022_2024/allshortreads/2022seq.Scaffold_9198_HRSCAF_10754.renamed.vcf \

 bgzip \
 -i -@ 10 \
  /scratch/rjp5nc/UK2022_2024/allshortreads/2022seq.Scaffold_9199_HRSCAF_10755.renamed.vcf \

 bgzip \
 -i -@ 10 \
  /scratch/rjp5nc/UK2022_2024/allshortreads/2022seq.Scaffold_9200_HRSCAF_10757.renamed.vcf \

 bgzip \
 -i -@ 10 \
  /scratch/rjp5nc/UK2022_2024/allshortreads/2022seq.Scaffold_9201_HRSCAF_10758.renamed.vcf \



 bcftools index \
  /scratch/rjp5nc/UK2022_2024/allshortreads/2022seq.Scaffold_1863_HRSCAF_2081.renamed.vcf.gz \

 bcftools index \
  /scratch/rjp5nc/UK2022_2024/allshortreads/2022seq.Scaffold_1931_HRSCAF_2197.renamed.vcf.gz \

 bcftools index \
  /scratch/rjp5nc/UK2022_2024/allshortreads/2022seq.Scaffold_2158_HRSCAF_2565.renamed.vcf.gz \

 bcftools index \
  /scratch/rjp5nc/UK2022_2024/allshortreads/2022seq.Scaffold_2217_HRSCAF_2652.renamed.vcf.gz \

 bcftools index \
  /scratch/rjp5nc/UK2022_2024/allshortreads/2022seq.Scaffold_2373_HRSCAF_2879.renamed.vcf.gz \

 bcftools index \
  /scratch/rjp5nc/UK2022_2024/allshortreads/2022seq.Scaffold_6786_HRSCAF_7541.renamed.vcf.gz \

 bcftools index \
  /scratch/rjp5nc/UK2022_2024/allshortreads/2022seq.Scaffold_7757_HRSCAF_8726.renamed.vcf.gz \

 bcftools index \
  /scratch/rjp5nc/UK2022_2024/allshortreads/2022seq.Scaffold_9197_HRSCAF_10753.renamed.vcf.gz \

 bcftools index \
  /scratch/rjp5nc/UK2022_2024/allshortreads/2022seq.Scaffold_9198_HRSCAF_10754.renamed.vcf.gz \

 bcftools index \
  /scratch/rjp5nc/UK2022_2024/allshortreads/2022seq.Scaffold_9199_HRSCAF_10755.renamed.vcf.gz \

 bcftools index \
  /scratch/rjp5nc/UK2022_2024/allshortreads/2022seq.Scaffold_9200_HRSCAF_10757.renamed.vcf.gz \

 bcftools index \
  /scratch/rjp5nc/UK2022_2024/allshortreads/2022seq.Scaffold_9201_HRSCAF_10758.renamed.vcf.gz \


bcftools concat \
--threads 10 \
-o  /scratch/rjp5nc/UK2022_2024/allshortreads/2022seq.concat.renamed.vcf.gz \
-O z \
/scratch/rjp5nc/UK2022_2024/allshortreads/2022seq.Scaffold_1863_HRSCAF_2081.renamed.vcf.gz \
/scratch/rjp5nc/UK2022_2024/allshortreads/2022seq.Scaffold_1931_HRSCAF_2197.renamed.vcf.gz \
/scratch/rjp5nc/UK2022_2024/allshortreads/2022seq.Scaffold_2217_HRSCAF_2652.renamed.vcf.gz \
/scratch/rjp5nc/UK2022_2024/allshortreads/2022seq.Scaffold_2158_HRSCAF_2565.renamed.vcf.gz \
/scratch/rjp5nc/UK2022_2024/allshortreads/2022seq.Scaffold_2373_HRSCAF_2879.renamed.vcf.gz \
/scratch/rjp5nc/UK2022_2024/allshortreads/2022seq.Scaffold_6786_HRSCAF_7541.renamed.vcf.gz \
/scratch/rjp5nc/UK2022_2024/allshortreads/2022seq.Scaffold_7757_HRSCAF_8726.renamed.vcf.gz \
/scratch/rjp5nc/UK2022_2024/allshortreads/2022seq.Scaffold_9197_HRSCAF_10753.renamed.vcf.gz \
/scratch/rjp5nc/UK2022_2024/allshortreads/2022seq.Scaffold_9198_HRSCAF_10754.renamed.vcf.gz \
/scratch/rjp5nc/UK2022_2024/allshortreads/2022seq.Scaffold_9199_HRSCAF_10755.renamed.vcf.gz \
/scratch/rjp5nc/UK2022_2024/allshortreads/2022seq.Scaffold_9200_HRSCAF_10757.renamed.vcf.gz \
/scratch/rjp5nc/UK2022_2024/allshortreads/2022seq.Scaffold_9201_HRSCAF_10758.renamed.vcf.gz
