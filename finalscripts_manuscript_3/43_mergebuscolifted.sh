#!/usr/bin/env bash

#SBATCH -J merge2 # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 24:00:00 ### 15 seconds
#SBATCH --mem 50G
#SBATCH -o /scratch/rjp5nc/erroroutputs/nFlo_1.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/erroroutputs/nFlo_1.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

module load bcftools 

cp /scratch/rjp5nc/UK2022_2024/buscoanalysis/filteredvcf/trimmed10bp_eupulex_vcf_old.vcf.gz_busco.vcf.gz /scratch/rjp5nc/UK2022_2024/buscoanalysis/filteredvcf_lifted/
cp /scratch/rjp5nc/UK2022_2024/buscoanalysis/filteredvcf/trimmed10bp_eupulex_vcf.vcf.gz_busco.vcf.gz /scratch/rjp5nc/UK2022_2024/buscoanalysis/filteredvcf_lifted/

tabix -p vcf /scratch/rjp5nc/UK2022_2024/buscoanalysis/filteredvcf_lifted/trimmed10bp_eupulex_vcf_old.vcf.gz_busco.vcf.gz
tabix -p vcf /scratch/rjp5nc/UK2022_2024/buscoanalysis/filteredvcf_lifted/trimmed10bp_eupulex_vcf.vcf.gz_busco.vcf.gz


cd /scratch/rjp5nc/UK2022_2024/buscoanalysis/filteredvcf_lifted/

bcftools merge -0 -Oz -o ../buscolifted_singlecopy.vcf.gz *.vcf.gz

#Assume homozygous ref - assume genotype is ref - use -0

tabix -p vcf ../buscolifted_singlecopy.vcf.gz

