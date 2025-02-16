#!/usr/bin/env bash
#
#SBATCH -J averagedepth # A single job name for the array
#SBATCH --ntasks-per-node=20 # ten cores
#SBATCH -N 1 # on one node
#SBATCH -t 16:00:00 ### 8 hours
#SBATCH --mem 100G
#SBATCH -o /scratch/rjp5nc/outputerrors/download_%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/outputerrors/download_%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

module load bcftools

# Define input and output files
VCF_FILE="/scratch/rjp5nc/UK2022_2024/allshortreads/chr/2022seq.concat.Removereps.renamed.vcf.gz"
OUTPUT_CSV="/scratch/rjp5nc/UK2022_2024/allshortreads/chr/read_depths.csv"

# Extract sample names for the header
SAMPLES=$(bcftools query -l "$VCF_FILE" | tr '\n' '\t')

# Add the header to the output file
echo -e "CHROM\tPOS\t${SAMPLES}" > "$OUTPUT_CSV"

# Extract read depth (DP) values and append to the file
bcftools query -f '%CHROM\t%POS[\t%DP]\n' -i 'DP>0' "$VCF_FILE" >> "$OUTPUT_CSV"