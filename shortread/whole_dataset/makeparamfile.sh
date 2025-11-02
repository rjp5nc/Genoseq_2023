#!/bin/bash

# Define the root directory where the search will start
root_dir="/scratch/rjp5nc/UK2022_2024/daphnia_phylo/gvcf/eupulex_chr"

# Output file where the results will be saved
output_file="/scratch/rjp5nc/UK2022_2024/daphnia_phylo/eu_pulexgvcfpaths.txt"

# Initialize the output file
echo -n > "$output_file"

# Loop through all .g.vcf.gz files in subdirectories
find "$root_dir" -type f -name "*.g.vcf.gz" | while read -r file; do
    # Extract the sample ID (assuming the format is sampleID.something.g.vcf.gz)
    sample_id=$(basename "$file" | cut -d'.' -f1)
    
    # Append the sample ID and file path to the output file
    echo -e "$sample_id\t$file" >> "$output_file"
done




cd /scratch/rjp5nc/UK2022_2024/daphnia_phylo/samplemapnames

while IFS= read -r line; do
  # Extract the basename of the directory (e.g., h2tg000008l)
  chr=$(echo "$line" | awk -F'/' '{for(i=1;i<=NF;i++) if($i ~ /^h2tg[0-9]+[lc]$/) print $i}')

  # Write the line to the corresponding .txt file
  echo "$line" >> "${chr}.txt"
done < "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/gvcf3/paths.txt"










cd /scratch/rjp5nc/UK2022_2024/daphnia_phylo/samplemapnames

while IFS= read -r line; do
  # Extract the chromosome identifier using grep-compatible regex
  chr=$(echo "$line" | grep -o 'JAACYE[0-9]\+\.[0-9]\+' | head -n1)

  # If an identifier was found, write the line to the corresponding file
  if [[ -n $chr ]]; then
    echo "$line" >> "${chr}.txt"
  else
    echo "Warning: No chromosome ID found in line: $line" >&2
  fi
done < "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/gvcf/usobtusa_chr/us_obtusapaths.txt"



cd /scratch/rjp5nc/UK2022_2024/daphnia_phylo/samplemapnames

input_file="/scratch/rjp5nc/UK2022_2024/daphnia_phylo/eu_pulexgvcfpaths.txt"

while IFS= read -r line; do
  # Extract scaffold/chromosome ID, assuming it appears in the path and has this pattern
  chr=$(echo "$line" | grep -oE 'Scaffold_[0-9]+_HRSCAF_[0-9]+' | head -n1)

  if [[ -n $chr ]]; then
    echo "$line" >> "${chr}.txt"
  else
    echo "Warning: No scaffold ID found in line: $line" >&2
  fi
done < "$input_file"






cat /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usambiguagvcfpaths.txt | head -n 5

cd /scratch/rjp5nc/UK2022_2024/daphnia_phylo/samplemapnames

input_file="/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usambiguagvcfpaths.txt"

while IFS= read -r line; do
  # Extract the chromosome ID (DAMBIG_3-style) from the path
  chr=$(echo "$line" | grep -oE 'DAMBIG_[0-9]+' | head -n1)

  if [[ -n $chr ]]; then
    echo "$line" >> "${chr}.txt"
  else
    echo "Warning: No chromosome ID found in line: $line" >&2
  fi
done < "$input_file"











cd /scratch/rjp5nc/UK2022_2024/daphnia_phylo/samplemapnames


input_file="/scratch/rjp5nc/UK2022_2024/daphnia_phylo/us_pulexgvcfpaths.txt"

while IFS= read -r line; do
  # Extract chromosome ID that looks like NC_ or JAACYE followed by numbers and a version
  chr=$(echo "$line" | grep -oE '(NC)_[0-9]+\.[0-9]+' | head -n1)

  # If chromosome ID is found, write the line to a corresponding file
  if [[ -n $chr ]]; then
    echo "$line" >> "${chr}.txt"
  else
    echo "Warning: No chromosome ID found in line: $line" >&2
  fi
done < "$input_file"
