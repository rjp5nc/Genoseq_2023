#!/usr/bin/env bash
#
#SBATCH -J CombineVCFs # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 5-00:00
#SBATCH --mem 100G
#SBATCH -o /scratch/rjp5nc/err/genotypegvcf.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/err/genotypegvcf.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab
#SBATCH --mail-type=END               # Send email at job completion
#SBATCH --mail-user=rjp5nc@virginia.edu    # Email address for notifications


# This script will merge all final VCFs together.

# Load modules
module load gatk/4.6.0.0

# Working folder is core folder where this pipeline is being run.

species=usobtusa_vcf

WORKING_FOLDER=/scratch/rjp5nc/UK2022_2024/daphnia_phylo/vcf/$species


# Intervals to analyze
# Create VCF combine list
cd $WORKING_FOLDER
ls *vcf.gz > $WORKING_FOLDER/${species}interval_paramList.list

# VCF list location
intervals=$WORKING_FOLDER/${species}interval_paramList.list

# Combined VCF name
output="/scratch/rjp5nc/UK2022_2024/daphnia_phylo/raw_vcf/raw_${species}.vcf.gz"

# Gatk parameters
JAVAMEM=40G
CPU=20

# Move to working directory
cd $WORKING_FOLDER

# Combine VCFs using MergeVcfs
gatk --java-options "-Xmx${JAVAMEM}" MergeVcfs \
-I $intervals \
-O $output

# Finish
echo "Complete" $(date)



#!/bin/bash

interval_list="$WORKING_FOLDER/${species}interval_paramList.list"  # Create this with the correct list of sample
reference_vcf="JAACYE010000216.1.1.51499.vcf.gz"

output_file="mismatched_vcfs.txt"

# Generate the expected sample list from the reference VCF
bcftools query -l "$reference_vcf" | sort > expected_samples.txt

# Initialize output
echo "VCFs with mismatched samples:" > "$output_file"

# Check each VCF in the list
while IFS= read -r vcf; do
    if [[ -f "$vcf" ]]; then
        echo "Checking $vcf..."
        bcftools query -l "$vcf" | sort > current_samples.txt
        if ! diff -q expected_samples.txt current_samples.txt > /dev/null; then
            echo "$vcf" >> "$output_file"
        fi
    else
        echo "Warning: File not found - $vcf"
    fi
done < "$interval_list"

# Clean up
rm -f expected_samples.txt current_samples.txt

echo "Done. See '$output_file' for mismatches."


