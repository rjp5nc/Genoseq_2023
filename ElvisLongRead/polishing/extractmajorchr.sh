#!/usr/bin/env bash
#
#SBATCH -J subset # A single job name for the array
#SBATCH --ntasks-per-node=10 # ten core
#SBATCH -N 1 # on one node
#SBATCH -t 0-10:00  ### 48 hours
#SBATCH --mem 100G
#SBATCH -o /scratch/rjp5nc/err/basicstats.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/err/basicstats.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

module load bcftools

cd /scratch/rjp5nc/UK2022_2024/daphnia_phylo/trimmed_10bp_repeatmasked_vcf/lifted_vcf/

bcftools view -r Scaffold_1863_HRSCAF_2081,Scaffold_1931_HRSCAF_2197,Scaffold_2217_HRSCAF_2652,Scaffold_2158_HRSCAF_2565,Scaffold_2373_HRSCAF_2879,Scaffold_6786_HRSCAF_7541,Scaffold_7757_HRSCAF_8726,Scaffold_9197_HRSCAF_10753,Scaffold_9198_HRSCAF_10754,Scaffold_9199_HRSCAF_10755,Scaffold_9200_HRSCAF_10757,Scaffold_9201_HRSCAF_10758 \
  -Oz -o lifted_12major.vcf.gz lifted_all.vcf.gz



