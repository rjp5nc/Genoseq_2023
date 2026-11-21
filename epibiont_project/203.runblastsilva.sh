#!/usr/bin/env bash
#SBATCH --job-name=blast_ssu_L6
#SBATCH --output=/scratch/rjp5nc/ssu_analysis/blast_test/blast_L6_%j.out
#SBATCH --error=/scratch/rjp5nc/ssu_analysis/blast_test/blast_L6_%j.err
#SBATCH --time=4:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --ntasks=1
#SBATCH --partition=standard
#SBATCH --account=berglandlab_standard
#SBATCH --array=0-1

set -euo pipefail

module load blast/2.17.0

BASE=/scratch/rjp5nc/ssu_analysis
DB="$BASE/databases/blast/silva_ssu"

LANES=(
    Gilmer5_A1_CKDL250003065-1A_22M5YKLT4_L3
    Gilmer5_A1_CKDL250003065-1A_22M5YKLT4_L4
)

LANE="${LANES[$SLURM_ARRAY_TASK_ID]}"
OUTDIR="$BASE/blast_test/$LANE"

blastn \
    -query "$OUTDIR/aligned_unique.fasta" \
    -db "$DB" \
    -task megablast \
    -num_threads "$SLURM_CPUS_PER_TASK" \
    -max_target_seqs 10 \
    -evalue 1e-10 \
    -outfmt "6 qseqid sseqid pident length qlen qcovhsp evalue bitscore stitle" \
    -out "$OUTDIR/aligned_all_hits.tsv"

echo "Finished $LANE"
wc -l "$OUTDIR/aligned_all_hits.tsv"
