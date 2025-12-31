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



# OUTDIR="/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/gatk_gvcf"
# mkdir -p "$OUTDIR"

# DIFFCSV="/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/mito_types_diff_df.csv"
# BAMDIR="/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/final_mitobam_rg3"

# awk -F',' 'NR>1 {gsub(/"/,"",$0); if ($4 > 10000) print $2}' "$DIFFCSV" | sort -u > "$OUTDIR/keep.samples"

# : > "$OUTDIR/bams.sitesgt10000.list"
# while read -r s; do
#   bam="${BAMDIR}/${s}finalmap_RG.bam"
#   [ -f "$bam" ] && echo "$bam" >> "$OUTDIR/bams.sitesgt10000.list"
# done < "$OUTDIR/keep.samples"




module load gatk samtools

OUTDIR="/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/gatk_gvcf"
BAMLIST="$OUTDIR/bams.sitesgt10000.list"
REF="/scratch/rjp5nc/Reference_genomes/mito_reference/usdobtusa_mito.fasta"
PLOIDY=2   # set to 1 if you want haploid mito

mkdir -p "$OUTDIR/gvcf" "$OUTDIR/logs" "$OUTDIR/tmp"

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











# module load gatk bcftools htslib samtools

# OUTDIR="/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/gatk_gvcf"
# REF="/scratch/rjp5nc/Reference_genomes/mito_reference/usdobtusa_mito.fasta"

# mkdir -p "$OUTDIR/tmp"

# # Ensure dict exists
# DICT="${REF%.fasta}.dict"
# [ -f "$DICT" ] || gatk CreateSequenceDictionary -R "$REF" -O "$DICT"

# # Build a file-of-filenames list for CombineGVCFs
# ls -1 "$OUTDIR/gvcf/"*.g.vcf.gz > "$OUTDIR/gvcfs.list"
# echo "gVCFs: $(wc -l < "$OUTDIR/gvcfs.list")"

# COMBINED="$OUTDIR/usdobtusa_mito_combined_from_bams.g.vcf.gz"
# JOINT="$OUTDIR/usdobtusa_mito_joint.vcf.gz"
# BI="$OUTDIR/usdobtusa_mito_biallelic.vcf.gz"
# CLEAN="$OUTDIR/usdobtusa.mito.biallelic.clean.vcf"

# # CombineGVCFs (fine for mito-sized data)
# gatk --java-options "-Xmx24g -Djava.io.tmpdir=$OUTDIR/tmp" CombineGVCFs \
#   -R "$REF" \
#   $(awk '{print "-V", $1}' "$OUTDIR/gvcfs.list") \
#   -O "$COMBINED"

# # Joint genotype
# gatk --java-options "-Xmx24g -Djava.io.tmpdir=$OUTDIR/tmp" GenotypeGVCFs \
#   -R "$REF" \
#   -V "$COMBINED" \
#   -O "$JOINT"

# # Biallelic SNPs
# bcftools view --threads "$SLURM_CPUS_PER_TASK" -m2 -M2 -v snps -Oz -o "$BI" "$JOINT"
# bcftools index -f "$BI"

# # Clean ALT="*"
# bcftools view -e 'ALT="*"' -m2 -M2 -Ov -o "$CLEAN" "$BI"

# echo "Final clean VCF -> $CLEAN"
# echo "Samples: $(bcftools query -l "$CLEAN" | wc -l)"
