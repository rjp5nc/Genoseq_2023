#!/usr/bin/env bash

#SBATCH -J BEAST # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 6-0:00:00 ### 15 seconds
#SBATCH --mem 60G
#SBATCH -o /scratch/rjp5nc/erroroutputs/beast.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/erroroutputs/beast.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab
#SBATCH --mail-type=END               # Send email at job completion
#SBATCH --mail-user=rjp5nc@virginia.edu    # Email address for notifications

# # Optional: create a new environment
# conda create -n snpeff_env openjdk=11 -y






# Load Modules
module load bcftools
module load gcc/9.2.0 htslib/1.10.2
module load tabix/0.2.6

conda activate snpeff_env

SNPEFF_HOME="/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/snpeff"
GENOME="US_DOBTUS"

MITO_FASTA="/scratch/rjp5nc/Reference_genomes/mito_reference/usdobtusa_mito.fasta"
MITOS_GFF="/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/snpeff/genes.gff"


bcftools norm -m -any -Oz -o /scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/snpeff/mito.split.vcf.gz /scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/usdobtusa_mito_allsites_all.haploid.vcf.gz
bcftools index /scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/snpeff/mito.split.vcf.gz

VCF_IN="/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/snpeff/mito.split.vcf.gz"

OUT_GFF="${SNPEFF_HOME}/genes.snpeff.gff3"
OUT_VCF="${SNPEFF_HOME}/output.ann.vcf"

# -------------------------
# ENV
# -------------------------
# (optional; comment out if you already activated it)
# conda activate snpeff_env

mkdir -p "${SNPEFF_HOME}/data/${GENOME}"

# -------------------------
# 1) Make snpEff-friendly GFF3
#    Convert protein-coding exons to CDS for mito genes.
#    (nad, cox, atp, cob/cytb). Keeps other features unchanged.
# -------------------------
awk -F'\t' 'BEGIN{OFS="\t"}
  /^#/ {print; next}
  NF < 9 {print; next}
  {
    feat=$3
    attr=$9
    parent=""

    # Extract Parent=... if present
    if (match(attr, /(^|;)Parent=([^;]+)/, m)) parent=m[2]

    # Convert exons that belong to protein-coding mito transcripts into CDS
    if (feat=="exon" && parent ~ /^transcript_(nad|cox|atp|cob|cytb)/) {
      $3="CDS"

      # phase must be 0/1/2 for CDS; if missing, set 0 (snpEff will still work)
      if ($8=="." || $8=="") $8="0"

      # ensure an ID= exists for CDS features
      if ($9 !~ /(^|;)ID=/) {
        newid="cds_" parent "_" $4 "_" $5
        $9=$9 ";ID=" newid
      }
    }

    print
  }
' "${MITOS_GFF}" > "${OUT_GFF}"

echo "Wrote: ${OUT_GFF}"
echo "First CDS lines:"
grep -n $'\tCDS\t' "${OUT_GFF}" | head

# -------------------------
# 2) Install genome files where snpEff expects them
# -------------------------
ln -sf "${MITO_FASTA}" "${SNPEFF_HOME}/data/${GENOME}/sequences.fa"
ln -sf "${OUT_GFF}"    "${SNPEFF_HOME}/data/${GENOME}/genes.gff"

echo "== Files in genome dir =="
ls -lah "${SNPEFF_HOME}/data/${GENOME}"

# -------------------------
# 3) Build snpEff.config starting from the DEFAULT config
#    (so codon tables are defined), then append our genome settings
#    IMPORTANT: use '=' (no spaces)
# -------------------------
DEFAULT_CFG=$(ls -1 "$CONDA_PREFIX"/share/snpeff-*/snpEff.config | head -n 1)
echo "Using default config: ${DEFAULT_CFG}"
cp -f "${DEFAULT_CFG}" "${SNPEFF_HOME}/snpEff.config"

cat >> "${SNPEFF_HOME}/snpEff.config" <<EOF

# ---- Custom genome ----
data_dir=${SNPEFF_HOME}/data
${GENOME}.genome=${GENOME}

# Invertebrate mitochondrial codon table
# (RHS must match codon table NAME in config, WITHOUT 'codon.' prefix)
${GENOME}.CM028013.1.codonTable=Invertebrate_Mitochondrial
EOF

echo "== Custom config lines =="
grep -n -E "data_dir=|${GENOME}\.genome=|CM028013\.1\.codonTable=" "${SNPEFF_HOME}/snpEff.config" | tail -n 20

# -------------------------
# 4) Build database
# -------------------------
snpEff build -c "${SNPEFF_HOME}/snpEff.config" -gff3 -noCheckCds -noCheckProtein "${GENOME}"

test -f "${SNPEFF_HOME}/data/${GENOME}/snpEffectPredictor.bin"
echo "OK: Built ${SNPEFF_HOME}/data/${GENOME}/snpEffectPredictor.bin"

# -------------------------
# 5) Annotate VCF
# -------------------------
snpEff -c "${SNPEFF_HOME}/snpEff.config" "${GENOME}" "${VCF_IN}" > "${OUT_VCF}"
echo "OK: Wrote ${OUT_VCF}"

# -------------------------
# 6) Quick verification: show missense/synonymous hits
# -------------------------
echo "== Example missense/synonymous hits =="
grep -v "^#" "${OUT_VCF}" | head -n 5000 | grep -E "missense_variant|synonymous_variant" | head -n 20 || true

# -------------------------
# 7) Count missense vs synonymous
#    (counts records that contain those strings anywhere in ANN)
# -------------------------
echo "== Counts =="
echo -n "missense_variant: "
grep -v "^#" "${OUT_VCF}" | grep -c "missense_variant" || true
echo -n "synonymous_variant: "
grep -v "^#" "${OUT_VCF}" | grep -c "synonymous_variant" || true



