#!/usr/bin/env bash
#
#SBATCH -J downloadSRA
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH -N 1
#SBATCH -t 0-48:00
#SBATCH --mem=20G
#SBATCH -o /scratch/rjp5nc/err/downloadsra.%A_%a.out
#SBATCH -e /scratch/rjp5nc/err/downloadsra.%A_%a.err
#SBATCH -p standard
#SBATCH --array=1-$(wc -l < /scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/sra_USobtusa.txt)
#SBATCH --account berglandlab

module load gcc/11.4.0 sratoolkit/3.1.1

# Paths
SRA_LIST="/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/sra_USobtusa.txt"
BASE_OUTDIR="/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/fasta"

# Get SRA accession for this array task
SRA=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SRA_LIST")

echo "Processing $SRA on task ${SLURM_ARRAY_TASK_ID}"

# Create output directory
OUTDIR="${BASE_OUTDIR}/${SRA}"
mkdir -p "$OUTDIR"
cd "$OUTDIR"

# Download
echo "Downloading $SRA..."
prefetch --max-size 200G "$SRA"

# Convert to FASTQ
echo "Converting $SRA to FASTQ..."
fasterq-dump "$SRA" \
  --split-files \
  --threads "${SLURM_CPUS_PER_TASK}" \
  --outdir "$OUTDIR" \
  --progress

# Compress
echo "Compressing FASTQ files..."
gzip -f "${SRA}"_*.fastq

echo "$SRA completed successfully."