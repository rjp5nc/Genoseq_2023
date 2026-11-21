#!/usr/bin/env bash
#
#SBATCH -J kraken_big
#SBATCH --ntasks-per-node=16
#SBATCH -N 1
#SBATCH -t 1:00:00
#SBATCH --mem 40G
#SBATCH -o /scratch/rjp5nc/err/kraken_big.%A_%a.out
#SBATCH -e /scratch/rjp5nc/err/kraken_big.%A_%a.err
#SBATCH -p standard
#SBATCH --account berglandlab_standard


set -euo pipefail

# ============================================================
# DIRECTORIES
# ============================================================


# # wget --timeout=30 --tries=1 \
#    -O SILVA_138_SSURef_NR99_tax_silva.fasta.gz \
#    https://ftp.arb-silva.de/release_138/Exports/SILVA_138_SSURef_NR99_tax_silva.fasta.gz
# --2026-11-20 21:11:24--  https://ftp.arb-silva.de/release_138/Exports/SILVA_138_SSURef_NR99_tax_silva.fasta.gz

#SBATCH --time=08:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4

set -euo pipefail

module load sortmerna/2.1

BASE=/scratch/rjp5nc/ssu_analysis
REF=$BASE/databases/silva_ssu.fasta
IDX=$BASE/databases/index/silva_ssu

mkdir -p "$BASE/databases/index"
mkdir -p "$BASE/databases/tmp"

indexdb_rna \
  --ref "${REF},${IDX}" \
  --tmpdir "$BASE/databases/tmp" \
  -m 16384 \
  -L 18 \
  -v