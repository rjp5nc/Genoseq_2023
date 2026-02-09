#!/usr/bin/env bash
#SBATCH -J DownloadMap    # Job name
#SBATCH --ntasks=1        # Single task per job
#SBATCH --cpus-per-task=10 # Number of CPU cores per task
#SBATCH -N 1              # Run on one node
#SBATCH -t 1:20:00        # 10 hours runtime
#SBATCH --mem=20G        # Memory per node
#SBATCH -o /scratch/rjp5nc/outputerrors/down.%A_%a.out  # Standard output
#SBATCH -e /scratch/rjp5nc/outputerrors/down.%A_%a.err  # Standard error
#SBATCH -p standard       # Partition
#SBATCH --account=berglandlab
#SBATCH --array=1-426%100   # Adjust based on the number of samples


#cat /scratch/rjp5nc/erroroutputs/beast.6896936_4294967294.out


module load gatk samtools

OUTDIR="/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/gatk_gvcf_hap"
BAMLIST="$OUTDIR/bams.sitesgt10000.list"
REF="/scratch/rjp5nc/Reference_genomes/mito_reference/usdobtusa_mito.fasta"
PLOIDY=1   

mkdir -p "$OUTDIR/gvcf" "$OUTDIR/logs" "$OUTDIR/tmp"



# DIFFCSV="/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/mito_types_diff_df.csv"
# BAMDIR="/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/final_mitobam_rg3"

# awk -F',' 'NR>1 {gsub(/"/,"",$0); if ($4 > 1000) print $2}' "$DIFFCSV" | sort -u > "$OUTDIR/keep.samples"

# : > "$OUTDIR/bams.sitesgt10000.list"
# while read -r s; do
#   bam="${BAMDIR}/${s}finalmap_RG.bam"
#   [ -f "$bam" ] && echo "$bam" >> "$OUTDIR/bams.sitesgt10000.list"
# done < "$OUTDIR/keep.samples"



# Grab this task's BAM
bam="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$BAMLIST")"
if [ -z "${bam:-}" ]; then
  echo "ERROR: No BAM at line ${SLURM_ARRAY_TASK_ID} in $BAMLIST" >&2
  exit 1
fi
if [ ! -f "$bam" ]; then
  echo "ERROR: BAM not found: $bam" >&2
  exit 1
fi

sample="$(basename "$bam")"
sample="${sample%finalmap_RG.bam}"   # adjust if suffix differs

# Ensure BAM index exists
[ -f "${bam}.bai" ] || samtools index "$bam"

# Ensure dict exists (GATK needs it)
DICT="${REF%.fasta}.dict"
[ -f "$DICT" ] || gatk CreateSequenceDictionary -R "$REF" -O "$DICT"

outg="$OUTDIR/gvcf/${sample}.g.vcf.gz"

gatk --java-options "-Xmx12g -Djava.io.tmpdir=$OUTDIR/tmp" HaplotypeCaller \
  -R "$REF" \
  -I "$bam" \
  -O "$outg" \
  -ERC GVCF \
  --sample-ploidy "$PLOIDY" \
  --native-pair-hmm-threads "$SLURM_CPUS_PER_TASK"

echo "Wrote $outg"










