#!/usr/bin/env bash
#SBATCH -J silva_index_v7
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH -N 1
#SBATCH -t 08:00:00
#SBATCH --mem=40G
#SBATCH -o /scratch/rjp5nc/err/silva_index_v7.%j.out
#SBATCH -e /scratch/rjp5nc/err/silva_index_v7.%j.err
#SBATCH -p standard
#SBATCH --account=berglandlab_standard

set -euo pipefail

# Initialize Conda for non-interactive SLURM jobs
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate sortmerna7

which sortmerna
sortmerna --version

BASE=/scratch/rjp5nc/ssu_analysis
REF=$BASE/databases/silva_ssu.fasta

# Separate from the old SortMeRNA 2.1 index
IDX=$BASE/databases/index_v7
WORK=$BASE/index_work_v7

mkdir -p "$IDX" "$WORK"

echo "Starting SILVA index: $(date)"

sortmerna \
    --ref "$REF" \
    --idx-dir "$IDX" \
    --workdir "$WORK" \
    --task 5 \
    --threads "${SLURM_CPUS_PER_TASK}" \
    -m 16384 \
    -L 18

echo "Indexing finished: $(date)"

du -sh "$IDX"
ls -lh "$IDX"

# # wget --timeout=30 --tries=1 \
#    -O SILVA_138_SSURef_NR99_tax_silva.fasta.gz \
#    https://ftp.arb-silva.de/release_138/Exports/SILVA_138_SSURef_NR99_tax_silva.fasta.gz
# --2026-11-20 21:11:24--  https://ftp.arb-silva.de/release_138/Exports/SILVA_138_SSURef_NR99_tax_silva.fasta.gz