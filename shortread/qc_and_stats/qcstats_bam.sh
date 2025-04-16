#!/usr/bin/env bash
#
#SBATCH -J qcstats # A single job name for the array
#SBATCH --ntasks-per-node=20 # ten cores
#SBATCH -N 1 # on one node
#SBATCH -t 16:00:00 ### 8 hours
#SBATCH --mem 100G
#SBATCH -o /scratch/rjp5nc/outputerrors/download_%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/outputerrors/download_%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

module load picard
module load samtools

cd /scratch/rjp5nc/UK2022_2024/mapped_bam

# Output CSV header
echo "Sample,Total_Reads,Mapped_Reads,Percent_Mapped,Duplicate_Rate" > bam_qc_summary.csv

# Loop through all BAM files recursively
find . -type f -name "*.bam" | while read bam; do
  sample=$(basename "$bam" .bam)

  # Total reads
  total_reads=$(samtools view -c "$bam")

  # Mapped reads (exclude unmapped, i.e., flag 0x4)
  mapped_reads=$(samtools view -c -F 4 "$bam")

  # Percent mapped
  if [ "$total_reads" -gt 0 ]; then
    percent_mapped=$(awk "BEGIN {printf \"%.2f\", ($mapped_reads / $total_reads) * 100}")
  else
    percent_mapped="0.00"
  fi

  # Duplicate rate via Picard (to temp file)
  tmp_metrics=$(mktemp)
  java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
    I="$bam" \
    O=/dev/null \
    M="$tmp_metrics" \
    REMOVE_DUPLICATES=false \
    ASSUME_SORTED=true \
    VALIDATION_STRINGENCY=SILENT \
    QUIET=true

  dup_line=$(grep -m1 'Unknown Library' "$tmp_metrics")
  dup_rate=$(echo "$dup_line" | awk '{if ($10 > 0) printf "%.4f", $9 / $10; else print "NA"}')

  rm "$tmp_metrics"

  # Append to CSV
  echo "$sample,$total_reads,$mapped_reads,$percent_mapped,$dup_rate" >> bam_qc_summary.csv
done