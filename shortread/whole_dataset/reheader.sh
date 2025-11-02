


#!/bin/bash

# Input file
input="/scratch/rjp5nc/UK2022_2024/mapped_bam/bam_files.txt"

# Output file path (change if needed)
output="/scratch/rjp5nc/UK2022_2024/mapped_bam/bam_files_with_names.txt"

# Create output folder if it doesn't exist
mkdir -p "$(dirname "$output")"

# Process each line
while IFS= read -r line; do
  # Extract filename
  fname=$(basename "$line")

  # Remove the suffix to get the sample name
  sample=$(echo "$fname" | sed 's/_finalmap\.bam$//')

  # Output original line + sample name
  echo -e "${line}\t${sample}"
done < "$input" > "$output"



rm *.vcf.gz.tbi


#!/bin/bash
# Path to the mapping file (full BAM path <TAB> sample name)
mapping_file="/scratch/rjp5nc/UK2022_2024/mapped_bam/bam_files_with_names.txt"

# Output directory
OUTDIR="/scratch/rjp5nc/UK2022_2024/renamed_vcfs"
mkdir -p "$OUTDIR"

module load bcftools

# Loop through each g.vcf.gz file
for file in /scratch/rjp5nc/UK2022_2024/daphnia_phylo/gvcf/usobtusa_chr/JAACYE010000084.1/*.g.vcf.gz; do
  # Extract current sample name (usually the full BAM path)
  current_sample=$(bcftools query -l "$file")

  # Find the matching line in the mapping file
  new_sample=$(awk -v curr="$current_sample" '$1 == curr {print $2}' "$mapping_file")

  if [[ -n "$new_sample" ]]; then
    base=$(basename "$file")
    out="$OUTDIR/$base"

    echo "Renaming sample in $file"
    echo -e "$current_sample\t$new_sample" | bcftools reheader -s /dev/stdin -o "$out" "$file"

    tabix -p vcf "$out"
  else
    echo "Warning: No mapping found for $current_sample in $file" >&2
  fi
done