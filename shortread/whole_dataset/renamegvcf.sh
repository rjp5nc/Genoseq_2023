#!/usr/bin/env bash
#
#SBATCH -J rename # A single job name for the array
#SBATCH --ntasks-per-node=10 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 12:00:00 # 8 hours
#SBATCH --mem 10G
#SBATCH -o /scratch/rjp5nc/err/gatk.chrom.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/err/gatk.chrom.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# This script will merge gVCFs into a unified database for genotype calling.
# This will be done using a per chromosome approach

target_dir="/scratch/rjp5nc/UK2022_2024/daphnia_phylo/gvcf"

# Find and rename files with carriage returns in their names
find "$target_dir" -type f | while IFS= read -r file; do
  if [[ "$file" == *$'\r'* ]]; then
    clean_file="$(echo "$file" | tr -d '\r')"
    echo "Renaming: $file -> $clean_file"
    mv "$file" "$clean_file"
  fi
done
