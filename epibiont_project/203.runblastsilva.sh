#!/usr/bin/env bash
#SBATCH --job-name=blast_ssu
#SBATCH --output=/scratch/rjp5nc/ssu_analysis/blast_test/blast_%j.out
#SBATCH --error=/scratch/rjp5nc/ssu_analysis/blast_test/blast_%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --ntasks=1
#SBATCH --partition=standard
#SBATCH --account=berglandlab_standard

set -euo pipefail

module load blast/2.17.0

BASE=/scratch/rjp5nc/ssu_analysis

echo "Started: $(date)"

blastn \
  -query "$BASE/blast_test/first_10000.fasta" \
  -db "$BASE/databases/blast/silva_ssu" \
  -task megablast \
  -num_threads "${SLURM_CPUS_PER_TASK}" \
  -max_target_seqs 10 \
  -evalue 1e-10 \
  -outfmt "6 qseqid sseqid pident length qlen qstart qend evalue bitscore stitle" \
  -out "$BASE/blast_test/first_10000_blast.tsv"

echo "Finished: $(date)"