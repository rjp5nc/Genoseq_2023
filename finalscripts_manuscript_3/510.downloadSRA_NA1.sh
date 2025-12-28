#!/usr/bin/env bash
#
#SBATCH -J downloadSRA # A single job name for the array
#SBATCH --ntasks-per-node=5 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-48:00  ### 48 hours
#SBATCH --mem 20G
#SBATCH -o /scratch/rjp5nc/downloadsra.out # Standard output
#SBATCH -e /scratch/rjp5nc/downloadsra.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

SRA_ACCESSION="SRR10160568"
OUTDIR="/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/fasta"

mkdir -p "$OUTDIR"
cd "$OUTDIR" || exit 1

echo "Downloading $SRA_ACCESSION..."
prefetch --max-size 200G "$SRA_ACCESSION"

echo "Converting to FASTQ..."
fasterq-dump "$SRA_ACCESSION" \
  --split-files \
  --threads 10 \
  --outdir "$OUTDIR" \
  --progress

echo "Compressing FASTQ files..."
gzip "${SRA_ACCESSION}"_*.fastq

echo "$SRA_ACCESSION completed."


