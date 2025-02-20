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

module load samtools

output_file="/scratch/rjp5nc/bamoutall.txt"
echo "Filename Average_Depth" > $output_file  # Header for the output file

for bam in /scratch/rjp5nc/UK2022_2024/allshortreads/sortedbamsdedup/*.bam; do 
    avg_depth=$(samtools depth -a "$bam" | awk '{sum+=$3} END {if (NR>0) print sum/NR; else print 0}')
    echo "$(basename "$bam") $avg_depth" >> $output_file
done