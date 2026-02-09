#!/usr/bin/env bash
#
#SBATCH -J downloadSRA # A single job name for the array
#SBATCH --ntasks-per-node=5 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-48:00  ### 48 hours
#SBATCH --mem 20G
#SBATCH -o /scratch/rjp5nc/err/downloadsra.out # Standard output
#SBATCH -e /scratch/rjp5nc/err/downloadsra.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

module load gcc/11.4.0 sratoolkit/3.1.1

SRA="SRR10160568"
OUTDIR="/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/fasta"
mkdir -p "$OUTDIR"
cd "$OUTDIR"

echo "Tools:"
which prefetch
which fasterq-dump

echo "Downloading $SRA..."
prefetch --max-size 200G "$SRA"

echo "Converting..."
fasterq-dump "$SRA" --split-files --threads 10 --outdir "$OUTDIR" --progress

echo "Compressing..."
gzip -f "${SRA}"_*.fastq

echo "Done."

