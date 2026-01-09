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


module load gatk samtools

META=/scratch/rjp5nc/UK2022_2024/touseforDBI_mito_fullref.csv
BAMDIR=/scratch/rjp5nc/UK2022_2024/final_mitobam_rg2

OUTDIR=/scratch/rjp5nc/UK2022_2024/euobtusa_mito/allsites_mito/gatk_gvcf_hap
REF=/scratch/rjp5nc/Reference_genomes/mito_reference/eudobtusa_mito_reverse.fasta
PLOIDY=1

mkdir -p "$OUTDIR/gvcf" "$OUTDIR/logs" "$OUTDIR/tmp"
cd "$OUTDIR" || exit 1

# ------------------------------------------------------------
# 1) Pull European Daphnia obtusa accessions from META
# ------------------------------------------------------------
awk -F',' 'NR>1 && $2=="Daphnia obtusa" && $3=="Europe" {print $1}' "$META" \
  | sort -u > EU_Daphnia_obtusa_accessions.txt

# ------------------------------------------------------------
# 2) Build BAM list for those accessions by searching BAMDIR
#    (matches anywhere in filename)
# ------------------------------------------------------------
# > bams.list
# while read -r acc; do
#   # grab first matching bam (if you expect exactly one)
#   bam=$(ls "$BAMDIR"/*"$acc"*".bam" 2>/dev/null | head -n 1)
#   if [[ -n "$bam" ]]; then
#     echo "$bam" >> bams.list
#   else
#     echo "MISSING_BAM $acc" >> missing_bams.txt
#   fi
# done < EU_Daphnia_obtusa_accessions.txt

# # If you want the job to fail when anything is missing:
# if [[ -f missing_bams.txt ]]; then
#   echo "Found missing BAMs (see $OUTDIR/missing_bams.txt). Failing." >&2
#   exit 2
# fi


# cd /scratch/rjp5nc/UK2022_2024/euobtusa_mito/allsites_mito/gatk_gvcf

# # remove empty/whitespace-only lines + deduplicate while preserving order
# awk 'NF' bams.list | awk '!seen[$0]++' > bams.list.clean
# mv bams.list.clean bams.list


# ------------------------------------------------------------
# 3) Array size helper (print this once, then set --array accordingly)
# ------------------------------------------------------------
# echo "N_BAMS=$(wc -l < /scratch/rjp5nc/UK2022_2024/euobtusa_mito/allsites_mito/gatk_gvcf/bams.list)" >&2

# ------------------------------------------------------------
# 4) Grab this task's BAM
# ------------------------------------------------------------




bam="$(sed -n "${SLURM_ARRAY_TASK_ID}p" bams.list)"
if [[ -z "${bam:-}" ]]; then
  echo "ERROR: No BAM at line ${SLURM_ARRAY_TASK_ID} in bams.list" >&2
  exit 1
fi
if [[ ! -f "$bam" ]]; then
  echo "ERROR: BAM not found: $bam" >&2
  exit 1
fi

# sample name
sample="$(basename "$bam")"
sample="${sample%.bam}"

# Ensure BAM index exists
[[ -f "${bam}.bai" ]] || samtools index "$bam"

# Ensure dict exists for REF
DICT="${REF%.*}.dict"
[[ -f "$DICT" ]] || gatk CreateSequenceDictionary -R "$REF" -O "$DICT"

outg="$OUTDIR/gvcf/${sample}.g.vcf.gz"

gatk --java-options "-Xmx12g -Djava.io.tmpdir=$OUTDIR/tmp" HaplotypeCaller \
  -R "$REF" \
  -I "$bam" \
  -O "$outg" \
  -ERC GVCF \
  --sample-ploidy "$PLOIDY" \
  --native-pair-hmm-threads "${SLURM_CPUS_PER_TASK:-1}"

echo "Wrote $outg"



OUTDIR="/scratch/rjp5nc/UK2022_2024/euobtusa_mito/allsites_mito/gatk_gvcf"
cd "$OUTDIR"
find gvcf -name "*.g.vcf.gz" | sort > gvcfs.list
REF="/scratch/rjp5nc/Reference_genomes/mito_reference/eudobtusa_mito_reverse.fasta"
gatk --java-options "-Xmx32g" CombineGVCFs \
  -R "$REF" \
  $(sed 's/^/-V /' gvcfs.list) \
  -O euobtusa_mito.all_samples.hap.g.vcf.gz
tabix -p vcf euobtusa_mito.all_samples.hap.g.vcf.gz
