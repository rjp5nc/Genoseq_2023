#!/usr/bin/env bash

#SBATCH -J rename # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 6-0:00:00 ### 15 seconds
#SBATCH --mem 60G
#SBATCH -o /scratch/rjp5nc/erroroutputs/beast.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/erroroutputs/beast.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab
#SBATCH --mail-type=END               # Send email at job completion
#SBATCH --mail-user=rjp5nc@virginia.edu    # Email address for notifications

module load bcftools

#bcftools annotate --rename-chrs rename_chroms.txt input.vcf.gz -Oz -o renamed.vcf.gz

cd /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/

bcftools annotate --rename-chrs contigmap.txt trimmed10bp_allsites_usobtusa.vcf.gz -Oz -o trimmed10bp_allsites_usobtusa_renamed.vcf.gz
bcftools index trimmed10bp_allsites_usobtusa_renamed.vcf.gz

