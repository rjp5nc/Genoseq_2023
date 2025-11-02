#!/bin/bash

module load samtools

#ref=dambigua_mito
ref=eudobtusa_mito_reverse
#ref=eudpulex_mito
#ref=kap4Dpulex_mito
#ref=usdobtusa_mito

samtools faidx /scratch/rjp5nc/Reference_genomes/mito_reference/$ref.fasta

# Input reference genome index file
input_fai="/scratch/rjp5nc/Reference_genomes/mito_reference/$ref.fasta.fai"

# Output directory for BED files
output_dir="/scratch/rjp5nc/Reference_genomes/mito_reference"

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Loop through each scaffold in the reference genome
while read -r line; do
    # Extract scaffold name and length
    scaffold_name=$(echo "$line" | awk '{print $1}')
    scaffold_length=$(echo "$line" | awk '{print $2}')

    # Create a BED file for this scaffold
    echo -e "${scaffold_name}\t0\t${scaffold_length}" > "${output_dir}/${ref}.bed"

done < "$input_fai"

echo "BED files created in $output_dir"
