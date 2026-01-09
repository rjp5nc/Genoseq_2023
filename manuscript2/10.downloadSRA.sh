#!/usr/bin/env bash
#
#SBATCH -J download # A single job name for the array
#SBATCH --ntasks-per-node=5 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-48:00  ### 48 hours
#SBATCH --mem 20G
#SBATCH -o /scratch/rjp5nc/downloadsra.out # Standard output
#SBATCH -e /scratch/rjp5nc/downloadsra.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab
#SBATCH --array=1-622%50

cd /scratch/rjp5nc/rawdata/SRRsamps

# Get the SRA ID for this array task
SRA=$(sed -n "${SLURM_ARRAY_TASK_ID}p" /scratch/rjp5nc/rawdata/sra_ids.txt)

echo "[$(date)] Processing $SRA on task ${SLURM_ARRAY_TASK_ID}"

# Create directory
mkdir -p "$SRA"

# Download
prefetch "$SRA"

# Convert to FASTQ
fasterq-dump \
  --split-files \
  --threads ${SLURM_CPUS_PER_TASK} \
  --outdir "$SRA" \
  "$SRA"

# Compress
gzip "$SRA"/*.fastq
rm "$SRA"/*.sra

echo "[$(date)] $SRA completed"
