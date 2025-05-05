#!/bin/bash

# Define the root directory where the search will start
root_dir="/scratch/rjp5nc/gvcftest"

# Output file where the results will be saved
output_file="/scratch/rjp5nc/gvcftest/params.txt"

# Initialize the output file
echo -n > "$output_file"

# Loop through all .g.vcf.gz files in subdirectories
find "$root_dir" -type f -name "*.g.vcf.gz" | while read -r file; do
    # Extract the sample ID (assuming the format is sampleID.something.g.vcf.gz)
    sample_id=$(basename "$file" | cut -d'.' -f1)
    
    # Append the sample ID and file path to the output file
    echo -e "$sample_id\t$file" >> "$output_file"
done