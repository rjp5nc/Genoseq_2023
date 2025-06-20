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

cd /scratch/rjp5nc/UK2022_2024/daphnia_phylo/trimmed_10bp_repeatmasked_vcf/lifted_vcf/lifted_notcomb

bcftools merge -0 -Oz -o ../lifted_all_missingasref.vcf.gz *.vcf.gz

#Assume homozygous ref - assume genotype is ref - use -0

tabix -p ../lifted_all_missingasref.vcf.gz

