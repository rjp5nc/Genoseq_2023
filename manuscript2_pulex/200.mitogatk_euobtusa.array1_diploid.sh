#!/usr/bin/env bash

#SBATCH -J windowed_het # A single job name for the array
#SBATCH --ntasks-per-node=10 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-3:00:00 ### 15 seconds
#SBATCH --mem 30G
#SBATCH -o /scratch/rjp5nc/erroroutputs/beast.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/erroroutputs/beast.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab
#SBATCH --mail-type=END               # Send email at job completion
#SBATCH --mail-user=rjp5nc@virginia.edu    # Email address for notifications
#SBATCH --array=1-43

module load gatk samtools htslib

META="/scratch/rjp5nc/UK2022_2024/touseforDBI_mito_fullref.csv"
BAMDIR="/scratch/rjp5nc/UK2022_2024/final_mitobam_rg2"

OUTDIR="/scratch/rjp5nc/UK2022_2024/euobtusa_mito/allsites_mito/gatk_gvcf"
REF="/scratch/rjp5nc/Reference_genomes/mito_reference/eudobtusa_mito_reverse.fasta"
PLOIDY=2

mkdir -p "$OUTDIR/gvcf" "$OUTDIR/tmp"
cd "$OUTDIR"

# ------------------------------------------------------------
# 0) One-time prep (safe to run every time)
#    - build accession list
#    - build bams.list from BAMDIR using those accessions
#    - clean bams.list (remove blanks, dedup)
#    - create ref dict if missing
# ------------------------------------------------------------

ACC_TXT="$OUTDIR/EU_Daphnia_obtusa_accessions.txt"
BAMLIST="$OUTDIR/bams.list"
MISSING="$OUTDIR/missing_bams.txt"

# awk -F',' 'NR>1 && $2=="Daphnia obtusa" && $3=="Europe" {print $1}' "$META" \
#   | sort -u > "$ACC_TXT"

# : > "$BAMLIST"
# : > "$MISSING"
# while read -r acc; do
#   bam=$(ls "$BAMDIR"/*"${acc}"*".bam" 2>/dev/null | head -n 1 || true)
#   if [[ -n "${bam:-}" ]]; then
#     echo "$bam" >> "$BAMLIST"
#   else
#     echo "$acc" >> "$MISSING"
#   fi
# done < "$ACC_TXT"

# # remove empty/whitespace-only lines + deduplicate while preserving order
# awk 'NF' "$BAMLIST" | awk '!seen[$0]++' > "$BAMLIST.clean"
# mv "$BAMLIST.clean" "$BAMLIST"

# If you want the array to HARD FAIL when any BAM is missing, uncomment:
# if [[ -s "$MISSING" ]]; then
#   echo "ERROR: Missing BAMs listed in $MISSING" >&2
#   exit 2
# fi

# Ensure dict exists for REF (GATK needs it)
DICT="${REF%.fasta}.dict"
[[ -f "$DICT" ]] || gatk CreateSequenceDictionary -R "$REF" -O "$DICT"

# ------------------------------------------------------------
# 1) Grab this task's BAM (array index)
# ------------------------------------------------------------
bam="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$BAMLIST")"
if [[ -z "${bam:-}" ]]; then
  echo "ERROR: No BAM at line ${SLURM_ARRAY_TASK_ID} in $BAMLIST" >&2
  echo "N_BAMS=$(wc -l < "$BAMLIST")" >&2
  exit 1
fi
if [[ ! -f "$bam" ]]; then
  echo "ERROR: BAM not found: $bam" >&2
  exit 1
fi

sample="$(basename "$bam")"
sample="${sample%.bam}"

# Ensure BAM index exists
[[ -f "${bam}.bai" ]] || samtools index "$bam"

outg="$OUTDIR/gvcf/${sample}.g.vcf.gz"
outg_tbi="${outg}.tbi"

# ------------------------------------------------------------
# 2) (Re)generate this sample's gVCF cleanly
#    Remove any partial outputs first to avoid bgzip/tabix weirdness
# ------------------------------------------------------------
rm -f "$outg" "$outg_tbi"

gatk --java-options "-Xmx12g -Djava.io.tmpdir=$OUTDIR/tmp" HaplotypeCaller \
  -R "$REF" \
  -I "$bam" \
  -O "$outg" \
  -ERC GVCF \
  --sample-ploidy "$PLOIDY" \
  --native-pair-hmm-threads "${SLURM_CPUS_PER_TASK:-1}"

echo "Wrote $outg"

# ------------------------------------------------------------
# 3) Combine step (only run once, after the array is finished)
#    This block is safe, but should NOT be run by every array task.
#    Gate it to only run on the LAST array index.
# ------------------------------------------------------------
N_BAMS="$(wc -l < "$BAMLIST")"

if [[ "${SLURM_ARRAY_TASK_ID}" -eq "${N_BAMS}" ]]; then
  module load bcftools htslib

  COMBINED="$OUTDIR/euobtusa_mito.all_samples.g.vcf.gz"
  COMBINED_TBI="${COMBINED}.tbi"
  GVLIST="$OUTDIR/gvcfs.list"

  # build list of gvcfs that actually exist and are non-empty
  find "$OUTDIR/gvcf" -name "*.g.vcf.gz" -size +0c | sort > "$GVLIST"

  # remove old combined outputs so re-runs are clean
  rm -f "$COMBINED" "$COMBINED_TBI"

  gatk --java-options "-Xmx32g -Djava.io.tmpdir=$OUTDIR/tmp" CombineGVCFs \
    -R "$REF" \
    $(sed 's/^/-V /' "$GVLIST") \
    -O "$COMBINED"

  tabix -f -p vcf "$COMBINED"
  echo "Wrote combined cohort gVCF: $COMBINED"
  echo "Indexed: $COMBINED.tbi"
else
  echo "Not combining on task ${SLURM_ARRAY_TASK_ID}/${N_BAMS} (combine runs only on last task)."
fi





