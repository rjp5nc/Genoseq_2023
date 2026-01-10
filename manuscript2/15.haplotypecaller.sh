#!/usr/bin/env bash
#
#SBATCH -J gatk_chrom # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-10:00:00 # 8 hours
#SBATCH --mem 25G
#SBATCH -o /scratch/rjp5nc/err/gatk.chrom.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/err/gatk.chrom.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab
#SBATCH --array=1-9999%150
#SBATCH --mail-type=END               # Send email at job completion
#SBATCH --mail-user=rjp5nc@virginia.edu    # Email address for notifications
#SBATCH --array=1-651%60   # Adjust the range based on the number of folders





module load gatk/4.6.0.0
module load samtools/1.17
module load gcc/14.2.0 htslib/1.17
module load tabix/0.2.6

THREADS="${SLURM_CPUS_PER_TASK:-1}"
echo "[$(date)] JOB=${SLURM_JOB_ID:-NA} TASK=${SLURM_ARRAY_TASK_ID:-NA} THREADS=$THREADS"
which gatk bgzip tabix samtools || true


# -----------------------
# CONFIG
# -----------------------


# OUT="/scratch/rjp5nc/rawdata/dobtusa_all_ids.txt"

# cat \
#   /scratch/rjp5nc/rawdata/sra_metadata_out/sra_ids_dobtusa.txt \
#   /scratch/rjp5nc/rawdata/mysamps_ids_dobtusa_europe.txt \
# | sed 's/\r//g' \
# | awk 'NF{print $1}' \
# | sort -u \
# > "$OUT"

# echo "Wrote: $OUT"
# wc -l "$OUT"
# head "$OUT"

IDFILE="/scratch/rjp5nc/rawdata/dobtusa_all_ids.txt"

# input BAMs from your previous pipeline
BAMDIR="/scratch/rjp5nc/UK2022_2024/redone_mito/euobtusa/sortedbamsdedup_obtusa"

# IMPORTANT: reference to call against (mito reference, since these are mito-mapped BAMs)
REF="/scratch/rjp5nc/Reference_genomes/mito_reference/eudobtusa_mito_reverse.fasta"

# output gvcfs
OUTROOT="/scratch/rjp5nc/UK2022_2024/redone_mito/euobtusa/gvcf_obtusa"
mkdir -p "$OUTROOT"

# -----------------------
# GET SRR FOR THIS TASK
# -----------------------
samp=$(awk -v n="${SLURM_ARRAY_TASK_ID}" 'NR==n{gsub(/\r/,""); print $1; exit}' "${IDFILE}")

if [[ -z "${samp}" ]]; then
  echo "ERROR: empty SRR at task ${SLURM_ARRAY_TASK_ID} from ${IDFILE}" >&2
  echo "Not a SRR"
fi
if [[ ! "${samp}" =~ ^SRR[0-9]+$ ]]; then
  echo "ERROR: samp is not SRR-like: [$samp]" >&2
  echo "Line ${SLURM_ARRAY_TASK_ID}:" >&2
  sed -n "${SLURM_ARRAY_TASK_ID}p" "${IDFILE}" | cat -A >&2
  echo "working on $samp" 
fi

INBAM="${BAMDIR}/${samp}.sort.dedup.bam"
if [[ ! -s "$INBAM" ]]; then
  echo "ERROR: input BAM missing: $INBAM" >&2
  ls -lh "$BAMDIR" | head -n 30 >&2 || true
fi

# -----------------------
# Ensure REF indexes exist
# -----------------------
if [[ ! -s "$REF" ]]; then
  echo "ERROR: REF missing: $REF" >&2
fi

if [[ ! -s "${REF}.fai" ]]; then
  echo "[$(date)] samtools faidx $REF"
  samtools faidx "$REF"
fi

DICT="${REF%.*}.dict"
# If REF ends with .fasta, this makes .dict with same base name.
if [[ ! -s "$DICT" ]]; then
  echo "[$(date)] gatk CreateSequenceDictionary -R $REF"
  gatk CreateSequenceDictionary -R "$REF"
fi

# -----------------------
# HaplotypeCaller -> gVCF
# -----------------------
OUTVCF="${OUTROOT}/${samp}.g.vcf"
OUTGZ="${OUTVCF}.gz"

echo "[$(date)] Calling gVCF"
echo "  SAMPLE: $samp"
echo "  INBAM:  $INBAM"
echo "  REF:    $REF"
echo "  OUT:    $OUTVCF"

gatk --java-options "-Xmx24g" HaplotypeCaller \
  -R "$REF" \
  -I "$INBAM" \
  -O "$OUTVCF" \
  -ERC GVCF \
  --native-pair-hmm-threads "$THREADS"

# -----------------------
# bgzip + tabix
# -----------------------
echo "[$(date)] Compressing + indexing"
bgzip -f "$OUTVCF"
tabix -f -p vcf "$OUTGZ"

# quick sanity
if [[ ! -s "$OUTGZ" || ! -s "${OUTGZ}.tbi" ]]; then
  echo "ERROR: gVCF output missing after compress/index" >&2
fi

echo "[$(date)] Done: $samp"
echo "  $OUTGZ"
