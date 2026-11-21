#!/usr/bin/env bash
#SBATCH --job-name=silva_blast_db
#SBATCH --output=/scratch/rjp5nc/ssu_analysis/blast_test/build_db_%j.out
#SBATCH --error=/scratch/rjp5nc/ssu_analysis/blast_test/build_db_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --partition=standard
#SBATCH --account=berglandlab_standard

set -euo pipefail

module load blast/2.17.0

BASE=/scratch/rjp5nc/ssu_analysis

mkdir -p "$BASE/databases/blast"

makeblastdb \
  -in "$BASE/databases/silva_ssu.fasta" \
  -dbtype nucl \
  -parse_seqids \
  -out "$BASE/databases/blast/silva_ssu"

echo "Database build completed: $(date)"