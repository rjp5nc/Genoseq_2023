#!/bin/bash

module load samtools

ref=assembly.hap2_onlydaps

samtools faidx /scratch/rjp5nc/Reference_genomes/post_kraken/$ref.fasta

# Input reference genome index file
input_fai="/scratch/rjp5nc/Reference_genomes/post_kraken/$ref.fasta.fai"

# Output directory for BED files
output_dir="/scratch/rjp5nc/Reference_genomes/$ref/scaffold_bed_files"

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Loop through each scaffold in the reference genome
while read -r line; do
    # Extract scaffold name and length
    scaffold_name=$(echo "$line" | awk '{print $1}')
    scaffold_length=$(echo "$line" | awk '{print $2}')

    # Create a BED file for this scaffold
    echo -e "${scaffold_name}\t0\t${scaffold_length}" > "${output_dir}/${scaffold_name}.bed"

done < "$input_fai"

echo "BED files created in $output_dir"











output_bed="/scratch/rjp5nc/Reference_genomes/${ref}/${ref}_all_scaffolds.bed"

# Create output directory if needed
mkdir -p "$(dirname "$output_bed")"

# Create a single BED file with all scaffolds
awk '{print $1"\t0\t"$2}' "$input_fai" > "$output_bed"

echo "Combined BED file created: $output_bed"