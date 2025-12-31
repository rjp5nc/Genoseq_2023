cd /scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/
module purge
module load gatk bcftools htslib samtools ruby

BAMDIR="/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/final_mitobam_rg3"
REF="/scratch/rjp5nc/Reference_genomes/mito_reference/usdobtusa_mito.fasta"
DIFFCSV="/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/mito_types_diff_df.csv"

OUTDIR="/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/gatk_gvcf"
mkdir -p "$OUTDIR"/{gvcf,logs,tmp}

# snapp_prep.rb (optional at end)
SNAPP_RB="/scratch/rjp5nc/snapp5/snapp_prep.rb"

# threshold
SITESTHRESH=10000

# -----------------------
# Reference index + dict (GATK needs dict)
# -----------------------
[ -f "${REF}.fai" ] || samtools faidx "$REF"
DICT="${REF%.fasta}.dict"
[ -f "$DICT" ] || gatk CreateSequenceDictionary -R "$REF" -O "$DICT"

# -----------------------
# 1) Sample list: compared_sites > 10000 (ignore diffs)
# -----------------------
KEEP_SAMPLES="$OUTDIR/keep.samples"
awk -F',' -v s="$SITESTHRESH" '
NR>1 { gsub(/"/,"",$0); if ($4 > s) print $2 }
' "$DIFFCSV" | sort -u > "$KEEP_SAMPLES"

echo "Samples passing compared_sites > ${SITESTHRESH}: $(wc -l < "$KEEP_SAMPLES")"

# -----------------------
# 2) Build BAM list (existing only)
# -----------------------
BAMLIST="$OUTDIR/bams.sitesgt${SITESTHRESH}.list"
: > "$BAMLIST"

while read -r s; do
  bam="${BAMDIR}/${s}finalmap_RG.bam"
  if [ -f "$bam" ]; then
    echo "$bam" >> "$BAMLIST"
  else
    echo "WARN missing BAM: $s ($bam)" >&2
  fi
done < "$KEEP_SAMPLES"

echo "BAMs found: $(wc -l < "$BAMLIST")"
head "$BAMLIST" || true

# -----------------------
# 3) Per-sample gVCF calling (HaplotypeCaller -ERC GVCF)
# -----------------------
# NOTE:
# - Assumes BAMs are coordinate-sorted and indexed (.bai). If not, index them.
#   If you want haploid mitochondria, set --sample-ploidy 1.
PLOIDY=1

while read -r bam; do
  sample="$(basename "$bam")"
  sample="${sample%finalmap_RG.bam}"   # yields sample name prefix
  outg="$OUTDIR/gvcf/${sample}.g.vcf.gz"

  # Index BAM if needed
  [ -f "${bam}.bai" ] || samtools index "$bam"

  echo "Calling gVCF for $sample"

  gatk --java-options "-Xmx8g -Djava.io.tmpdir=$OUTDIR/tmp" HaplotypeCaller \
    -R "$REF" \
    -I "$bam" \
    -O "$outg" \
    -ERC GVCF \
    --sample-ploidy "$PLOIDY" \
    > "$OUTDIR/logs/${sample}.haplotypecaller.out" \
    2> "$OUTDIR/logs/${sample}.haplotypecaller.err"

done < "$BAMLIST"

# Verify gVCFs
ls -1 "$OUTDIR/gvcf/"*.g.vcf.gz | wc -l

# -----------------------
# 4) Combine gVCFs (simple CombineGVCFs approach)
# -----------------------
# For large N, GenomicsDBImport is better, but CombineGVCFs is fine for mito-sized data.
COMBINED_GVCF="$OUTDIR/usdobtusa_mito_combined_from_bams.g.vcf.gz"

# Build a list of "-V file" args for GATK
VCF_ARGS=()
while read -r f; do
  VCF_ARGS+=("-V" "$f")
done < <(ls -1 "$OUTDIR/gvcf/"*.g.vcf.gz)

gatk --java-options "-Xmx16g -Djava.io.tmpdir=$OUTDIR/tmp" CombineGVCFs \
  -R "$REF" \
  "${VCF_ARGS[@]}" \
  -O "$COMBINED_GVCF"

# -----------------------
# 5) Joint genotyping
# -----------------------
JOINT="$OUTDIR/usdobtusa_mito_joint.vcf.gz"
gatk --java-options "-Xmx16g -Djava.io.tmpdir=$OUTDIR/tmp" GenotypeGVCFs \
  -R "$REF" \
  -V "$COMBINED_GVCF" \
  -O "$JOINT"

# -----------------------
# 6) Biallelic SNPs + clean ALT="*"
# -----------------------
BI="$OUTDIR/usdobtusa_mito_biallelic.vcf.gz"
bcftools view --threads 15 -m2 -M2 -v snps -Oz -o "$BI" "$JOINT"
bcftools index -f "$BI"

CLEAN="$OUTDIR/usdobtusa.mito.biallelic.clean.vcf"
bcftools view -e 'ALT="*"' -m2 -M2 -Ov -o "$CLEAN" "$BI"

echo "Final clean VCF for SNAPP -> $CLEAN"
echo "Samples in clean VCF: $(bcftools query -l "$CLEAN" | wc -l)"

# -----------------------
# 7) (Optional) SNAPP prep files + XML
# -----------------------
SAMPLES="$OUTDIR/samples_mito_USobtusa_snapp_clean.txt"
MAP2COL="$OUTDIR/samples_mito_USobtusa_snapp_clean_2col.txt"
CONSTRAINTS="$OUTDIR/snapp_constraints.txt"
XML="$OUTDIR/snapp.mito.xml"

bcftools query -l "$CLEAN" > "$SAMPLES"

echo -n "monophyletic NA " > "$CONSTRAINTS"
sed 's/$/_clone/' "$SAMPLES" | paste -sd, - >> "$CONSTRAINTS"

awk '{print $1 "_clone", $1}' "$SAMPLES" > "$MAP2COL"

ruby "$SNAPP_RB" \
  -v "$CLEAN" \
  -t "$MAP2COL" \
  -c "$CONSTRAINTS" \
  -x "$XML" \
  -o "$OUTDIR/snapp.mito" \
  -m 1000 \
  -l 10000

echo "SNAPP XML -> $XML"