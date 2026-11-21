
#!/usr/bin/env bash
#SBATCH --job-name=sortmerna7_L3L4
#SBATCH --output=/scratch/rjp5nc/ssu_analysis/logs/sortmerna7_L3L4_%A_%a.out
#SBATCH --error=/scratch/rjp5nc/ssu_analysis/logs/sortmerna7_L3L4_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --mem=40G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=standard
#SBATCH --account=berglandlab_standard
#SBATCH --array=1-576%30

set -euo pipefail

source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate sortmerna7

BASE=/scratch/rjp5nc/ssu_analysis
REF="$BASE/databases/silva_ssu.fasta"
IDX="$BASE/databases/index_v7"
SAMPLES="$BASE/samples.txt"
OUTBASE="$BASE/results_v7"

TASK_ID=${SLURM_ARRAY_TASK_ID}
SAMPLE_DIR=$(sed -n "${TASK_ID}p" "$SAMPLES")

if [[ -z "$SAMPLE_DIR" || ! -d "$SAMPLE_DIR" ]]; then
    echo "ERROR: Invalid sample directory: $SAMPLE_DIR"
    exit 1
fi

SAMPLE=$(basename "$SAMPLE_DIR")
OUT="$OUTBASE/$SAMPLE"

mkdir -p "$OUT"

echo "Sample: $SAMPLE"
echo "Started: $(date)"

# Select L3 and L4 sequencing lanes only
shopt -s nullglob
R1_FILES=(
    "$SAMPLE_DIR"/*_L3_1.fq.gz
    "$SAMPLE_DIR"/*_L4_1.fq.gz
)

if (( ${#R1_FILES[@]} == 0 )); then
    echo "No L3 or L4 lanes found for $SAMPLE. Skipping."
    exit 0
fi

for R1 in "${R1_FILES[@]}"; do

    R2="${R1%_1.fq.gz}_2.fq.gz"

    if [[ ! -f "$R2" ]]; then
        echo "ERROR: Missing mate: $R2"
        exit 1
    fi

    LANE=$(basename "${R1%_1.fq.gz}")
    WORK="$OUT/work_$LANE"

    echo "Processing lane: $LANE"
    echo "R1: $R1"
    echo "R2: $R2"

    mkdir -p "$WORK"

    sortmerna \
        --ref "$REF" \
        --reads "$R1" \
        --reads "$R2" \
        --idx-dir "$IDX" \
        --workdir "$WORK" \
        --task 4 \
        --fastx \
        --paired_in \
        --threads "${SLURM_CPUS_PER_TASK}"

    echo "Completed lane: $LANE"

done

echo "Completed sample: $SAMPLE"
echo "Finished: $(date)"
