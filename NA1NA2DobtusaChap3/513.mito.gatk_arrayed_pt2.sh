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


module load gatk bcftools htslib samtools

OUTDIR="/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/gatk_gvcf"
REF="/scratch/rjp5nc/Reference_genomes/mito_reference/usdobtusa_mito.fasta"



mkdir -p "$OUTDIR/tmp"

cd $OUTDIR

# Ensure dict exists
DICT="${REF%.fasta}.dict"
[ -f "$DICT" ] || gatk CreateSequenceDictionary -R "$REF" -O "$DICT"

# Build a file-of-filenames list for CombineGVCFs
ls -1 "$OUTDIR/gvcf/"*.g.vcf.gz > "$OUTDIR/gvcfs.list"
echo "gVCFs: $(wc -l < "$OUTDIR/gvcfs.list")"

COMBINED="$OUTDIR/usdobtusa_mito_combined_from_bams.g.vcf.gz"
JOINT="$OUTDIR/usdobtusa_mito_joint.vcf.gz"
FULL="$OUTDIR/usdobtusa_mito_joint_Full.vcf.gz"
BI="$OUTDIR/usdobtusa_mito_biallelic.vcf.gz"
CLEAN="$OUTDIR/usdobtusa.mito.biallelic.clean.vcf"
BIFULL="$OUTDIR/usdobtusa_mito_biallelic_FULL.vcf.gz"

# CombineGVCFs (fine for mito-sized data)
gatk --java-options "-Xmx24g -Djava.io.tmpdir=$OUTDIR/tmp" CombineGVCFs \
  -R "$REF" \
  $(awk '{print "-V", $1}' "$OUTDIR/gvcfs.list") \
  -O "$COMBINED"

# Joint genotype
gatk --java-options "-Xmx24g -Djava.io.tmpdir=$OUTDIR/tmp" GenotypeGVCFs \
  -R "$REF" \
  -V "$COMBINED" \
  -O "$JOINT"

# Joint genotype
gatk --java-options "-Xmx24g -Djava.io.tmpdir=$OUTDIR/tmp" GenotypeGVCFs \
  -R "$REF" \
  -V "$COMBINED" \
  --include-non-variant-sites \
  -O "$FULL"





# Biallelic SNPs
bcftools view --threads "$SLURM_CPUS_PER_TASK" -m2 -M2 -v snps -Oz -o "$BI" "$JOINT"
bcftools index -f "$BI"

bcftools view --threads "$SLURM_CPUS_PER_TASK" \
  -m1 -M2 \
  -Oz -o "$BIFULL" "$FULL"

bcftools index -f "$BIFULL"

zcat $BIFULL | wc -l

IN="$OUTDIR/usdobtusa_mito_biallelic.vcf.gz"
OUT="$OUTDIR/usdobtusa_mito_biallelic.diploidGT.vcf.gz"

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

BI="$OUTDIR/usdobtusa_mito_biallelic.diploidGT.vcf.gz"

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


#Haploid
#
# WARNING: The maximum number of SNPs has been set to 10000, which is greater
#     than the number of bi-allelic SNPs with sufficient information (371) for SNAPP.
# WARNING: Excluded 1691 sites with only missing data in one or more species.

# INFO: Retained 371 bi-allelic sites.