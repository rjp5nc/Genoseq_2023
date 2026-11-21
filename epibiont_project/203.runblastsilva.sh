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

set -euo pipefail

module load blast/2.17.0

BASE=/scratch/rjp5nc/ssu_analysis
SAMPLE=Gilmer5_A1
LANE=Gilmer5_A1_CKDL240041264-1A_22TG7VLT3_L6

OUTDIR="$BASE/blast_test/${LANE}"
mkdir -p "$OUTDIR"

ALIGNED="$BASE/results_v7/$SAMPLE/work_$LANE/out/aligned.fq.gz"

FASTA="$OUTDIR/aligned_unique.fasta"
BLAST="$OUTDIR/aligned_blast.tsv"
FILTERED="$OUTDIR/aligned_filtered.tsv"

echo "Starting: $(date)"
echo "Input: $ALIGNED"

# Step 1: Convert ALL aligned FASTQ records to FASTA.
# Preserve mate identity using /1 and /2.
gzip -cd "$ALIGNED" | \
awk '
  NR % 4 == 1 {
    split(substr($0, 2), h, " ")
    mate = substr(h[2], 1, 1)

    if (mate != "1" && mate != "2") {
      print "ERROR: Unexpected mate identifier: " $0 > "/dev/stderr"
      exit 1
    }

    print ">" h[1] "/" mate
  }

  NR % 4 == 2 {
    print
  }
' > "$FASTA"

# Step 2: Check total read count and unique identifiers
TOTAL=$(grep -c '^>' "$FASTA")
UNIQUE=$(grep '^>' "$FASTA" | sort -u | wc -l)

echo "Total reads: $TOTAL"
echo "Unique identifiers: $UNIQUE"

if [[ "$TOTAL" -ne "$UNIQUE" ]]; then
  echo "ERROR: Duplicate FASTA identifiers detected."
  exit 1
fi

# Step 3: BLAST ALL reads against SILVA
blastn \
    -query "$FASTA" \
    -db "$BASE/databases/blast/silva_ssu" \
    -task megablast \
    -num_threads "${SLURM_CPUS_PER_TASK:-8}" \
    -max_target_seqs 10 \
    -evalue 1e-10 \
    -outfmt "6 qseqid sseqid pident length qlen qcovhsp evalue bitscore stitle" \
    -out "$OUTDIR/aligned_blast_relaxed.tsv"

# Filter using relaxed thresholds
awk -F '\t' '
    $3 >= 90 && $6 >= 60
' "$OUTDIR/aligned_blast_relaxed.tsv" \
    > "$OUTDIR/aligned_filtered_relaxed.tsv"

# Step 5: Summary
echo "--------------------------------"
echo "BLAST summary"
echo "--------------------------------"

echo "Reads tested: $TOTAL"

echo -n "Reads with any BLAST hit: "
cut -f1 "$BLAST" | sort -u | wc -l

echo -n "Reads passing 95% identity / 80% coverage: "
cut -f1 "$FILTERED" | sort -u | wc -l

echo "Finished: $(date)"
echo "Filtered output: $FILTERED"
