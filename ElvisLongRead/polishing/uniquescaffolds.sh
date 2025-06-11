#!/usr/bin/env bash
#
#SBATCH -J subset # A single job name for the array
#SBATCH --ntasks-per-node=1 # ten core
#SBATCH -N 1 # on one node
#SBATCH -t 2-10:00  ### 48 hours
#SBATCH --mem 20G
#SBATCH -o /scratch/rjp5nc/err/basicstats.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/err/basicstats.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

module load bcftools

cd /scratch/rjp5nc/UK2022_2024/daphnia_phylo/trimmed_10bp_repeatmasked_vcf/lifted_vcf/
 
bcftools query -f '%CHROM\n' lifted_12major.vcf.gz | sort -u > unique_scaffolds.txt