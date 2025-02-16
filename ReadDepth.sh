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
VCF_FILE="/scratch/rjp5nc/UK2022_2024/allshortreads/chr/2022seq.concat.Removereps.renamed.vcf.gz"  # Replace with your actual VCF file path
OUTPUT_CSV="average_read_depth.csv"
# Extract read depth (DP) for each sample and compute the average
echo "Sample,Average_Read_Depth" > "$OUTPUT_CSV"
bcftools query -f '%CHROM\t%POS[\t%DP]\n' -i 'DP>0' /scratch/rjp5nc/UK2022_2024/allshortreads/chr/2022seq.concat.Removereps.renamed.vcf.gz --threads 10 > /scratch/rjp5nc/UK2022_2024/allshortreads/chr/read_depths.csv
