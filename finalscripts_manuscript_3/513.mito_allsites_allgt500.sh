#!/usr/bin/env bash
#SBATCH -J DownloadMap    # Job name
#SBATCH --ntasks=1        # Single task per job
#SBATCH --cpus-per-task=10 # Number of CPU cores per task
#SBATCH -N 1              # Run on one node
#SBATCH -t 1-20:00        # 10 hours runtime
#SBATCH --mem=50G        # Memory per node
#SBATCH -o /scratch/rjp5nc/outputerrors/down.%A_%a.out  # Standard output
#SBATCH -e /scratch/rjp5nc/outputerrors/down.%A_%a.err  # Standard error
#SBATCH -p standard       # Partition
#SBATCH --account=berglandlab



cd /scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/
module load bcftools samtools

# -----------------------
# Inputs
# -----------------------
BAMDIR="/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/final_mitobam_rg3"
REF="/scratch/rjp5nc/Reference_genomes/mito_reference/usdobtusa_mito.fasta"

DP="/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/usdobtusa_mito_avg_DP_per_sample.txt"
DIFFCSV="/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/mito_types_diff_df.csv"

OUTDIR="/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/vcf_split_by_diffs"
mkdir -p "$OUTDIR"

VCF_OUT="$OUTDIR/usdobtusa_mito_allsites_all.diploid.dpgt30.diffsgt500.vcf.gz"
BAMLIST="$OUTDIR/bams.dpgt30.diffsgt500.txt"

# -----------------------
# Reference index
# -----------------------
[ -f "${REF}.fai" ] || samtools faidx "$REF"

# -----------------------
# Build sample list: (DP>30) ∩ (diffs<500)
# -----------------------
awk '$2 > 30 {print $1}' "$DP" | sort -u > "$OUTDIR/dp_gt30.samples"
awk -F',' 'NR>1 {gsub(/"/,"",$0); if ($3 > 600 && $4 > 10000) print $2}' "$DIFFCSV" \
  | sort -u > "$OUTDIR/diffs_lt500.samples"
  comm -12 "$OUTDIR/dp_gt30.samples" "$OUTDIR/diffs_gt500.samples" > "$OUTDIR/keep.samples"

echo "Keeping samples: $(wc -l < "$OUTDIR/keep.samples")"

# -----------------------
# Turn sample list into BAM list (FULL PATHS), only if BAM exists
# -----------------------
: > "$BAMLIST"
while read -r s; do
  bam="${BAMDIR}/${s}finalmap_RG.bam"
  if [ -f "$bam" ]; then
    echo "$bam" >> "$BAMLIST"
  else
    echo "WARN missing BAM for sample: $s ($bam)" >&2
  fi
done < "$OUTDIR/keep.samples"

echo "BAMs found: $(wc -l < "$BAMLIST")"
head "$BAMLIST"

# -----------------------
# Call diploid VCF using only those BAMs
# -----------------------
bcftools mpileup \
  -f "$REF" \
  -q 20 -Q 20 \
  -a FORMAT/DP,FORMAT/AD \
  -Ou \
  -b "$BAMLIST" \
| bcftools call -m -A --ploidy 2 -Oz -o "$VCF_OUT"

bcftools index -f "$VCF_OUT"
echo "Wrote -> $VCF_OUT"
echo "VCF sample count: $(bcftools query -l "$VCF_OUT" | wc -l)"
