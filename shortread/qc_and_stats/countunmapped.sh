#!/usr/bin/env bash
#
#SBATCH -J countunmapped # A single job name for the array
#SBATCH --ntasks-per-node=20 # ten core
#SBATCH -N 1 # on one node
#SBATCH -t 0-12:00  ### 24 hours
#SBATCH --mem 100G
#SBATCH -o /scratch/rjp5nc/countunmapped.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/countunmapped.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

module load samtools

# Output file
output_file="/scratch/rjp5nc/UK2022_2024/allshortreads/unmappedinsortedbam.csv"

wd="/scratch/rjp5nc/UK2022_2024/allshortreads/sortedbams"
# Output CSV file

# Write header to CSV file
echo "BAM_File,Unmapped_Reads" > "$output_file"

# Loop through all BAM files in the current directory
for bam in *.bam; do
    # Count unmapped reads
    count=$(samtools view -c -f 4 "$bam")

    # Append results to CSV file
    echo "$bam,$count" >> "$output_file"
done

echo "Results saved in $output_file"