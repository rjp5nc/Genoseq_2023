#!/usr/bin/env bash

#SBATCH -J downloadSRA # A single job name for the array
#SBATCH --ntasks-per-node=10 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 2-0:00:00 ### 15 seconds
#SBATCH --mem 20G
#SBATCH -o /scratch/rjp5nc/err/downloadsra.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/err/downloadsra.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab
#SBATCH --mail-type=END               # Send email at job completion
#SBATCH --mail-user=rjp5nc@virginia.edu    # Email address for notifications


###sbatch --array=3-286%40 510.downloadSRA_other_obtusa.sh

module load gcc/11.4.0 sratoolkit/3.1.1

# Paths
SRA_LIST="/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/sra_USobtusa.txt"
BASE_OUTDIR="/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/fasta"

# Get SRA accession for this array task
SRA=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SRA_LIST")

echo "Processing $SRA on task ${SLURM_ARRAY_TASK_ID}"

# Create output directory
OUTDIR="${BASE_OUTDIR}"
mkdir -p "$OUTDIR"
cd "$OUTDIR"

# Download
echo "Downloading $SRA..."
prefetch --max-size 200G "$SRA"

# Convert to FASTQ
echo "Converting $SRA to FASTQ..."
fasterq-dump "${OUTDIR}/${SRA}/${SRA}.sra" \
  --split-files \
  --threads "${SLURM_CPUS_PER_TASK}" \
  --outdir "${OUTDIR}" \
  --progress

# Compress
echo "Compressing FASTQ files..."
gzip -f "${SRA}"_*.fastq

echo "$SRA completed successfully."