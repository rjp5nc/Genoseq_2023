#!/bin/bash

# Define input directory
INPUT_DIR="/scratch/rjp5nc/HMW/HMWDNAElvis3/blast/trimmedfasta"

# Count the number of FASTA files
NUM_FILES=$(ls "$INPUT_DIR"/*.fasta | wc -l)

# Submit the SLURM array job
sbatch --array=1-"$NUM_FILES" /home/rjp5nc/Genoseq_2023/ElvisLongRead/shblast.sh

echo "Submitted job array with $NUM_FILES tasks."





#after finished

cat /scratch/rjp5nc/HMW/HMWDNAElvis3/blast_results2/*_final_results.txt > /scratch/rjp5nc/HMW/HMWDNAElvis3/all_final_results.txt
