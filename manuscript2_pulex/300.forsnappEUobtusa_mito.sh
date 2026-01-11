


module load gatk bcftools htslib samtools vcftools ruby

OUTDIR="/scratch/rjp5nc/UK2022_2024/mitogvcf/gvcf/eudobtusa_mito_reverse_out"
REF="/scratch/rjp5nc/Reference_genomes/mito_reference/eudobtusa_mito_reverse.fasta"

# This is your *input* cohort gVCF (do not delete)
COMBINED="/scratch/rjp5nc/UK2022_2024/euobtusa_mito/allsites_mito/gatk_gvcf/euobtusa_mito.all_samples.g.vcf.gz"

mkdir -p "$OUTDIR/tmp"
cd "$OUTDIR"

# ----------------------------
# Ensure reference dict exists
# ----------------------------
DICT="${REF%.fasta}.dict"
if [[ ! -f "$DICT" ]]; then
  gatk CreateSequenceDictionary -R "$REF" -O "$DICT"
fi

# ----------------------------
# Define outputs (diploid pipeline products)
# ----------------------------
JOINT="$OUTDIR/eudobtusa_mito_joint.vcf.gz"
BI="$OUTDIR/eudobtusa_mito_biallelic.vcf.gz"
DIP="$OUTDIR/eudobtusa_mito_biallelic.diploidGT.vcf.gz"
CLEAN="$OUTDIR/eudobtusa.mito.biallelic.clean.vcf"

SAMPLES="$OUTDIR/samples_mito_USobtusa_snapp_clean.txt"
CONSTRAINTS="$OUTDIR/snapp_constraints.txt"
MAP2COL="$OUTDIR/samples_mito_USobtusa_snapp_clean_2col.txt"

SNAPP_XML="$OUTDIR/snapp.mito_dip.xml"
SNAPP_PREFIX="$OUTDIR/snapp.mito"   # snapp_prep writes multiple files with this prefix

# ----------------------------
# CLEANUP: remove diploid intermediates + downstream products
# ----------------------------
rm -f \
  "$JOINT" "$JOINT.tbi" "$JOINT.csi" \
  "$BI" "$BI.tbi" "$BI.csi" \
  "$DIP" "$DIP.tbi" "$DIP.csi" \
  "$CLEAN" \
  "$SAMPLES" "$CONSTRAINTS" "$MAP2COL" \
  "$SNAPP_XML"

# remove SNAPP outputs safely (only those created by this pipeline)
rm -f "${SNAPP_PREFIX}"*  || true

# vcftools missingness output (miss.imiss, miss.lmiss, miss.log)
rm -f miss.imiss miss.lmiss miss.log

# tmp can accumulate junk; optional wipe (uncomment if you want)
# rm -rf "$OUTDIR/tmp" && mkdir -p "$OUTDIR/tmp"

# ----------------------------
# 1) Joint genotype (diploid genotyping from combined gVCF)
# ----------------------------
gatk --java-options "-Xmx24g -Djava.io.tmpdir=$OUTDIR/tmp" GenotypeGVCFs \
  -R "$REF" \
  -V "$COMBINED" \
  -O "$JOINT"

tabix -f -p vcf "$JOINT"

# ----------------------------
# 2) Keep only biallelic SNPs
# ----------------------------
THREADS="${SLURM_CPUS_PER_TASK:-8}"
bcftools view --threads "$THREADS" -m2 -M2 -v snps -Oz -o "$BI" "$JOINT"
bcftools index -f "$BI"

# ----------------------------
# 3) Convert haploid GTs (0 or 1) to diploid (0/0 or 1/1)
#    ONLY when GT is a single allele "0" or "1" or "."
# ----------------------------
zcat "$BI" | awk 'BEGIN{OFS="\t"}
  /^#/ {print; next}
  {
    for(i=10;i<=NF;i++){
      n=split($i,a,":")
      # a[1] is GT
      if(a[1] ~ /^[0-9.]$/){
        a[1]=a[1] "/" a[1]
        $i=a[1]
        for(j=2;j<=n;j++) $i=$i ":" a[j]
      }
    }
    print
  }' | bgzip -c > "$DIP"

tabix -f -p vcf "$DIP"

# quick spot check
bcftools query -s SRR14370492 -f '%POS\t[%GT\t%DP\t%AD]\n' "$DIP" | head -n 30

# ----------------------------
# 4) Clean ALT="*" and write a plain VCF (SNAPP prep wants plain VCF often)
# ----------------------------
bcftools view -e 'ALT="*"' -m2 -M2 -Ov -o "$CLEAN" "$DIP"

echo "Final clean VCF -> $CLEAN"
echo "Samples: $(bcftools query -l "$CLEAN" | wc -l)"

# ----------------------------
# 5) Missingness report
# ----------------------------
vcftools --vcf "$CLEAN" --missing-indv --out miss
head miss.imiss

# ----------------------------
# 6) SNAPP prep
# ----------------------------
bcftools query -l "$CLEAN" > "$SAMPLES"

echo -n "monophyletic NA " > "$CONSTRAINTS"
sed 's/$/_clone/' "$SAMPLES" | paste -sd, - >> "$CONSTRAINTS"

awk '{print $1 "_clone", $1}' "$SAMPLES" > "$MAP2COL"

ruby /scratch/rjp5nc/snapp5/snapp_prep.rb \
  -v "$CLEAN" \
  -t "$MAP2COL" \
  -c "$CONSTRAINTS" \
  -x "$SNAPP_XML" \
  -o "$SNAPP_PREFIX" \
  -m 10000 \
  -l 10000