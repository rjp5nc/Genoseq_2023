#!/bin/bash

# Define input directory
INPUT_DIR="/scratch/rjp5nc/HMW/HMWDNAElvis3/blast/trimmed_fasta"

# Count the number of FASTA files
NUM_FILES=$(ls "$INPUT_DIR"/*.fasta | wc -l)

# Submit the SLURM array job
sbatch --array=1-"$NUM_FILES" blast_worker.sh

echo "Submitted job array with $NUM_FILES tasks."