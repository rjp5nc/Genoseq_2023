

module load gatk bcftools htslib samtools

OUTDIR="/scratch/rjp5nc/UK2022_2024/mitogvcf/gvcf/eudobtusa_mito_reverse_out"
REF="/scratch/rjp5nc/Reference_genomes/mito_reference/eudobtusa_mito_reverse.fasta"

mkdir -p "$OUTDIR/tmp"

cd $OUTDIR

# Ensure dict exists
DICT="${REF%.fasta}.dict"
[ -f "$DICT" ] || gatk CreateSequenceDictionary -R "$REF" -O "$DICT"

# Build a file-of-filenames list for CombineGVCFs
ls -1 "/scratch/rjp5nc/UK2022_2024/mitogvcf/gvcf/eudobtusa_mito_reverse/"*.g.vcf.gz > "$OUTDIR/gvcfs.list"
echo "gVCFs: $(wc -l < "$OUTDIR/gvcfs.list")"

COMBINED="$/scratch/rjp5nc/UK2022_2024/euobtusa_mito/allsites_mito/gatk_gvcf/euobtusa_mito.all_samples.g.vcf.gz"
JOINT="$OUTDIR/eudobtusa_mito_joint.vcf.gz"
BI="$OUTDIR/eudobtusa_mito_biallelic.vcf.gz"
CLEAN="$OUTDIR/eudobtusa.mito.biallelic.clean.vcf"

# Joint genotype
gatk --java-options "-Xmx24g -Djava.io.tmpdir=$OUTDIR/tmp" GenotypeGVCFs \
  -R "$REF" \
  -V "$COMBINED" \
  -O "$JOINT"

# Biallelic SNPs
bcftools view --threads "$SLURM_CPUS_PER_TASK" -m2 -M2 -v snps -Oz -o "$BI" "$JOINT"
bcftools index -f "$BI"



IN="$OUTDIR/eudobtusa_mito_biallelic.vcf.gz"
OUT="$OUTDIR/eudobtusa_mito_biallelic.diploidGT.vcf.gz"

zcat "$IN" | awk 'BEGIN{OFS="\t"}
  /^#/ {print; next}
  {
    for(i=10;i<=NF;i++){
      n=split($i,a,":")
      if(a[1] ~ /^[0-9.]$/){
        a[1]=a[1] "/" a[1]
        $i=a[1]
        for(j=2;j<=n;j++) $i=$i ":" a[j]
      }
    }
    print
  }' | bgzip -c > "$OUT"

tabix -p vcf "$OUT"

BI="$OUTDIR/eudobtusa_mito_biallelic.diploidGT.vcf.gz"

# Clean ALT="*"
bcftools view -e 'ALT="*"' -m2 -M2 -Ov -o "$CLEAN" "$BI"






echo "Final clean VCF -> $CLEAN"
echo "Samples: $(bcftools query -l "$CLEAN" | wc -l)"




module load bcftools htslib ruby

# ---- sample list from clean VCF ----
SAMPLES="$OUTDIR/samples_mito_USobtusa_snapp_clean.txt"
bcftools query -l "$CLEAN" > "$SAMPLES"

# ---- constraints + 2-col mapping ----
CONSTRAINTS="$OUTDIR/snapp_constraints.txt"
MAP2COL="$OUTDIR/samples_mito_USobtusa_snapp_clean_2col.txt"

echo -n "monophyletic NA " > "$CONSTRAINTS"
sed 's/$/_clone/' "$SAMPLES" | paste -sd, - >> "$CONSTRAINTS"

awk '{print $1 "_clone", $1}' "$SAMPLES" > "$MAP2COL"

# ---- SNAPP prep ----
ruby /scratch/rjp5nc/snapp5/snapp_prep.rb \
  -v "$CLEAN" \
  -t "$MAP2COL" \
  -c "$CONSTRAINTS" \
  -x "$OUTDIR/snapp.mito_dip.xml" \
  -o "$OUTDIR/snapp.mito" \
  -m 10000 \
  -l 10000
