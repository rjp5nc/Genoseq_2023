#!/usr/bin/env bash
#
#SBATCH -J pixy # A single job name for the array
#SBATCH --ntasks-per-node=10 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-92:00  ### 48 hours
#SBATCH --mem 80G
#SBATCH -o /scratch/rjp5nc/err/pixy.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/err/pixy.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab
module load bcftools

cd /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv
mkdir -p chr

# Get unique chromosomes
chroms=$(bcftools query -f '%CHROM\n' trimmed10bp_allsites_usobtusa.filtered_bgz.vcf.gz | sort -u)

# Loop over each chromosome
for chr in $chroms; do
    echo "Processing $chr"

    # Make sure the VCF is indexed first
    bcftools index trimmed10bp_allsites_usobtusa.filtered_bgz.vcf.gz

    # Extract chromosome-specific VCF
    bcftools view -r $chr -Oz -o chr/"${chr}_invariant.vcf.gz" trimmed10bp_allsites_usobtusa.filtered_bgz.vcf.gz

    # Index the new chromosome VCF
    bcftools index -t chr/"${chr}_invariant.vcf.gz"
done
