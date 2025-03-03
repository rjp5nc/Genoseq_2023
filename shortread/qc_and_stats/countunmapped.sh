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

# sbatch ~/Genoseq_2023/shortread/qc_and_stats/countunmapped.sh


module load samtools

# Output file
output_file="/scratch/rjp5nc/UK2022_2024/allshortreads/unmappedinsortedbam.csv"

wd="/scratch/rjp5nc/UK2022_2024/allshortreads/sortedbams"
cd ${wd}

echo "BAM_File,Mapped,Unmapped_Reads" > "$output_file"

# Output CSV file

for bam in *.bam; do
    # Count mapped reads
    mapped=$(samtools view -c -F 4 "$bam")

    # Count unmapped reads
    unmapped=$(samtools view -c -f 4 "$bam")

    # Count total reads
    total=$(samtools view -c "$bam")

    # Append results to CSV file
    echo "$bam,$mapped,$unmapped,$total" >> "$output_file"
done
