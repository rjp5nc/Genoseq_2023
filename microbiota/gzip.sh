#!/usr/bin/env bash
#
#SBATCH -J gzip # A single job name for the array
#SBATCH --ntasks-per-node=5 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-3:00  ### 48 hours
#SBATCH --mem 100G
#SBATCH -o /scratch/rjp5nc/downloadsra.out # Standard output
#SBATCH -e /scratch/rjp5nc/downloadsra.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab


SRA_ACCESSION="SRR14476950"

# Create output directory
OUTPUT_DIR="/scratch/rjp5nc/microbiota"

# Compress FASTQ files
echo "Compressing FASTQ files..."
gzip "$OUTPUT_DIR/${SRA_ACCESSION}_1.fastq"
gzip "$OUTPUT_DIR/${SRA_ACCESSION}_2.fastq"
