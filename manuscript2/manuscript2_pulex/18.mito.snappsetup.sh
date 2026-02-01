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

module load bcftools htslib samtools ruby

# -----------------------
# CONFIG (EU obtusa all-sites cohort VCF)
# -----------------------
IN_ALLSITES="/scratch/rjp5nc/UK2022_2024/redone_mito/euobtusa/cohort_gendb/obtusa.mito.ALLSITES.rmMissGT0.3.vcf.gz"

OUTDIR="/scratch/rjp5nc/UK2022_2024/redone_mito/euobtusa/snapp_from_allsites"
mkdir -p "$OUTDIR"
cd "$OUTDIR"

# -----------------------
# sanity
# -----------------------
if [[ ! -s "$IN_ALLSITES" ]]; then
  echo "ERROR: missing input VCF: $IN_ALLSITES" >&2
  exit 1
fi

# index if needed (either .tbi or .csi is fine)
if [[ ! -s "${IN_ALLSITES}.tbi" && ! -s "${IN_ALLSITES}.csi" ]]; then
  echo "Indexing input VCF..."
  bcftools index -f "$IN_ALLSITES"
fi

echo "Input all-sites VCF: $IN_ALLSITES"
echo "Samples: $(bcftools query -l "$IN_ALLSITES" | wc -l)"

# -----------------------
# 1) Keep only biallelic SNPs (variant sites only)
#    NOTE: this drops invariant sites, which is what SNAPP wants.
# -----------------------
BI="${OUTDIR}/euobtusa_mito_biallelic.vcf.gz"
bcftools view \
  --threads "${SLURM_CPUS_PER_TASK:-1}" \
  -m2 -M2 -v snps \
  -Oz -o "$BI" \
  "$IN_ALLSITES"
bcftools index -f "$BI"

# -----------------------
# 2) Convert haploid GT (0 or 1) to diploid (0/0 or 1/1)
#    (only touches GT field; keeps the rest)
# -----------------------
DIP="${OUTDIR}/euobtusa_mito_biallelic.diploidGT.vcf.gz"

zcat "$BI" | awk 'BEGIN{OFS="\t"}
  /^#/ {print; next}
  {
    for(i=10;i<=NF;i++){
      n=split($i,a,":")
      # a[1] is GT; handle 0,1,., and already-diploid cases
      if(a[1] ~ /^[01]$/){
        a[1]=a[1] "/" a[1]
      }
      $i=a[1]
      for(j=2;j<=n;j++) $i=$i ":" a[j]
    }
    print
  }' | bgzip -c > "$DIP"

tabix -f -p vcf "$DIP"

# -----------------------
# 3) Remove ALT="*" and keep biallelic SNPs (again, just to be safe)
# -----------------------
CLEAN="${OUTDIR}/euobtusa_mito_biallelic.diploidGT.clean.vcf"
bcftools view \
  -e 'ALT="*"' \
  -m2 -M2 -v snps \
  -Ov -o "$CLEAN" \
  "$DIP"

echo "Final clean VCF -> $CLEAN"
echo "Samples: $(bcftools query -l "$CLEAN" | wc -l)"

# -----------------------
# 4) SNAPP prep
# -----------------------
SAMPLES="${OUTDIR}/samples_mito_EUobtusa_snapp_clean.txt"
bcftools query -l "$CLEAN" > "$SAMPLES"

CONSTRAINTS="${OUTDIR}/snapp_constraints.txt"
MAP2COL="${OUTDIR}/samples_mito_EUobtusa_snapp_clean_2col.txt"

echo -n "monophyletic NA " > "$CONSTRAINTS"
sed 's/$/_clone/' "$SAMPLES" | paste -sd, - >> "$CONSTRAINTS"

awk '{print $1 "_clone", $1}' "$SAMPLES" > "$MAP2COL"




# ----------------------------
# (4) Run snapp_prep.rb -> SNAPP XML
# ----------------------------
XML_OUT="$OUTDIR/snapp.mito.euobtusa.xml"
PREFIX_OUT="$OUTDIR/snapp.mito.euobtusa"

ruby /scratch/rjp5nc/snapp5/snapp_prep.rb \
  -v "$CLEAN_VCF" \
  -t "$MAP_2COL" \
  -c "$CONSTR" \
  -x "$XML_OUT" \
  -o "$PREFIX_OUT" \
  -m 10000 \
  -l 10000
