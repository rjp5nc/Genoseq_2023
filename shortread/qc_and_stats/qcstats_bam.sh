#!/usr/bin/env bash
#
#SBATCH -J qcstats # A single job name for the array
#SBATCH --ntasks-per-node=5 # ten cores
#SBATCH -N 1 # on one node
#SBATCH -t 8:00:00 ### 8 hours
#SBATCH --mem 40G
#SBATCH -o /scratch/rjp5nc/outputerrors/download_%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/outputerrors/download_%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab
#SBATCH --array=1-519%40   # Adjust based on the number of samples

module load samtools
module load picard

cd /scratch/rjp5nc/UK2022_2024/mapped_bam

# Get BAM file for this task
BAM=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" bam_files.txt)
SAMPLE=$(basename "$BAM" .bam)

# QC Output directory
OUTDIR="qc_results"
mkdir -p "$OUTDIR"

# Total and mapped reads
TOTAL_READS=$(samtools view -c "$BAM")
MAPPED_READS=$(samtools view -c -F 4 "$BAM")
PCT_MAPPED=$(awk "BEGIN {printf \"%.2f\", ($MAPPED_READS / $TOTAL_READS) * 100}")

# Duplicate rate
TMP_METRICS=$(mktemp)
java -jar $EBROOTPICARD/picard.jar MarkDuplicates I="$BAM" O=/dev/null M="$TMP_METRICS" REMOVE_DUPLICATES=false ASSUME_SORTED=true QUIET=true VALIDATION_STRINGENCY=SILENT
DUP_LINE=$(grep -m1 'Unknown Library' "$TMP_METRICS")
DUP_RATE=$(echo "$DUP_LINE" | awk '{if ($10 > 0) printf "%.4f", $9 / $10; else print "NA"}')
rm "$TMP_METRICS"

# Average coverage (requires BAM index)
COVERAGE=$(samtools depth "$BAM" | awk '{sum += $3} END {if (NR > 0) printf "%.2f", sum / NR; else print "0.00"}')

# Output to per-sample CSV
echo "$SAMPLE,$TOTAL_READS,$MAPPED_READS,$PCT_MAPPED,$DUP_RATE,$COVERAGE" > "$OUTDIR/${SAMPLE}_qc.csv"


echo "Sample,Total_Reads,Mapped_Reads,Percent_Mapped,Duplicate_Rate,Average_Coverage" > bam_qc_summary.csv
cat qc_results/*_qc.csv >> bam_qc_summary.csv