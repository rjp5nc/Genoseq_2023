#!/usr/bin/env bash
#
#SBATCH -J trim # A single job name for the array
#SBATCH --ntasks-per-node=10 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-92:00  ### 48 hours
#SBATCH --mem 50G
#SBATCH -o /scratch/rjp5nc/err/trim.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/err/trim.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

module load bcftools


VCF="/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/trimmed10bp_allsites_Repeatmasked_usobtusa.filtered_bgz.vcf.gz"
SAMPLES="/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/heterozygotegilmer/onlygilmer.txt"
OUT="/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/trimmed10bp_filtered_Gilmer.vcf.gz"

# filter VCF by samples and bgzip output
bcftools view -S "$SAMPLES" -Oz -o "$OUT" "$VCF"

# index the new VCF
bcftools index "$OUT"
