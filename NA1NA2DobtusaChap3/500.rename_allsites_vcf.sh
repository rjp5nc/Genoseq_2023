#!/usr/bin/env bash

#SBATCH -J rename # A single job name for the array
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

module load bcftools

#bcftools annotate --rename-chrs rename_chroms.txt input.vcf.gz -Oz -o renamed.vcf.gz

# cd /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/



#cp trimmed10bp_allsites_usobtusa.vcf.gz trimmed10bp_allsites_usobtusa2.vcf.gz

# bcftools annotate --rename-chrs contigmap.txt trimmed10bp_allsites_usobtusa.vcf.gz -Oz -o trimmed10bp_allsites_usobtusa_renamed.vcf.gz
# bcftools index trimmed10bp_allsites_usobtusa_renamed.vcf.gz

# # awk -v MAP="contigmap.txt" '
# #   BEGIN{
# #     # read map; split on ANY whitespace
# #     while ((getline line < MAP) > 0) {
# #       if (line ~ /^[[:space:]]*$/) continue
# #       n = split(line, a, /[[:space:]]+/)
# #       if (n >= 2) m[a[1]] = a[2]
# #     }
# #     close(MAP)
# #   }
# #   # For FASTA
# #   /^>/{
# #     hdr = substr($0,2)                 # drop ">"
# #     id  = hdr
# #     if (match(hdr, /^[^[:space:]]+/))  # first token = true ID
# #       id = substr(hdr, RSTART, RLENGTH)
# #     rest = substr(hdr, length(id)+1)   # preserve description
# #     if (id in m) id = m[id]
# #     print ">" id rest
# #     next
# #   }
# #   { print }                            # sequence lines
# # ' /scratch/rjp5nc/Reference_genomes/post_kraken/US_obtusa_onlydaps.fa > US_obtusa_onlydaps.renamed.fa
# # cat US_obtusa_onlydaps.renamed.fa | head -n 5


# # awk '/^>/{sub(/_LG[0-9]+/,""); print; next} {print}' US_obtusa_onlydaps.renamed.fa > US_obtusa_onlydaps.nolgs.fa
# # samtools faidx US_obtusa_onlydaps.nolgs.fa


# # conda create -n gfftools -c bioconda -c conda-forge gffread samtools bcftools bedtools
# # conda activate gfftools

# awk 'BEGIN{OFS=FS="\t"} 
#      !/^#/ {
#         sub(/_LG[0-9]+$/, "", $1)
#      }
#      {print}
# ' Daphnia_obtusa_FS6_genome.gtf > Daphnia_obtusa_FS6_genome.renamed.gtf




#conda create -n snpeff -c bioconda -c conda-forge snpeff openjdk -y

cd /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/

conda activate snpeff
which snpEff



# ====== SETTINGS ======
GENOME_ID="US_DOBTUS"
REF_FASTA="US_obtusa_onlydaps.renamed.fa"
ANNOT_GTF="Daphnia_obtusa_FS6_genome.renamed.gtf"
INPUT_VCF="trimmed10bp_allsites_usobtusa_renamed.vcf.gz"
OUTPUT_VCF="trimmed10bp_allsites_usobtusa_renamed_annotated.vcf.gz"


GENOME_ID="US_DOBTUS"

REF_FASTA="US_obtusa_onlydaps.renamed.fa"
ANNOT_GTF="Daphnia_obtusa_FS6_genome.renamed.gtf"

INPUT_VCF="trimmed10bp_allsites_usobtusa_renamed.vcf.gz"
OUTPUT_VCF="trimmed10bp_allsites_usobtusa_renamed_annotated.vcf.gz"

# Locate snpEff install under this conda env
SNPEFF_HOME="$(find "${CONDA_PREFIX}/share" -maxdepth 1 -type d -name 'snpeff*' | head -n1)"
SNPEFF_CONFIG="${SNPEFF_HOME}/snpEff.config"
GENOME_DIR="${SNPEFF_HOME}/data/${GENOME_ID}"

echo "GENOME_ID   = $GENOME_ID"
echo "REF_FASTA   = $REF_FASTA"
echo "ANNOT_GTF   = $ANNOT_GTF"
echo "INPUT_VCF   = $INPUT_VCF"
echo "OUTPUT_VCF  = $OUTPUT_VCF"
echo "SNPEFF_HOME = $SNPEFF_HOME"
echo "SNPEFF_CONF = $SNPEFF_CONFIG"
echo "GENOME_DIR  = $GENOME_DIR"
echo

# ---- 2. Sanity checks on input files / config ----
[[ -f "$REF_FASTA"  ]] || { echo "ERROR: REF_FASTA not found: $REF_FASTA" >&2; exit 1; }
[[ -f "$ANNOT_GTF"  ]] || { echo "ERROR: ANNOT_GTF not found: $ANNOT_GTF" >&2; exit 1; }
[[ -f "$INPUT_VCF"  ]] || { echo "ERROR: INPUT_VCF not found: $INPUT_VCF" >&2; exit 1; }
[[ -f "$SNPEFF_CONFIG" ]] || { echo "ERROR: snpEff.config not found at $SNPEFF_CONFIG" >&2; exit 1; }

# ---- 3. Prepare genome directory (clean + copy GTF/FASTA) ----
echo "Setting up SnpEff data directory..."
rm -rf "$GENOME_DIR"
mkdir -p "$GENOME_DIR"

cp "$REF_FASTA"  "$GENOME_DIR/sequences.fa"
cp "$ANNOT_GTF"  "$GENOME_DIR/genes.gtf"

# ---- 4. Ensure genome is listed in snpEff.config ----
if ! grep -q "^${GENOME_ID}\.genome" "$SNPEFF_CONFIG"; then
  echo "Adding ${GENOME_ID} to snpEff.config"
  echo "${GENOME_ID}.genome : ${GENOME_ID}" >> "$SNPEFF_CONFIG"
else
  echo "Genome ${GENOME_ID} already present in snpEff.config"
fi

# ---- 5. Build SnpEff database from GTF ----
echo
echo "Building SnpEff database (skipping CDS/protein checks)..."

pushd "$SNPEFF_HOME" >/dev/null

# IMPORTANT: disable CDS and protein checks or DB will NOT be saved
java -Xmx4g -jar snpEff.jar build -gtf22 -v -noCheckCds -noCheckProtein "$GENOME_ID"

popd >/dev/null
echo "Build step finished."
echo

# Quick check that the DB exists now
if [[ ! -f "$GENOME_DIR/snpEffectPredictor.bin" ]]; then
  echo "ERROR: snpEffectPredictor.bin not found in $GENOME_DIR" >&2
  echo "Build failed; check the build log above." >&2
  exit 1
fi

# ---- 6. Annotate VCF ----
echo "Annotating VCF with local database (no download)..."

# Remove any old output/index
rm -f "$OUTPUT_VCF" "$OUTPUT_VCF.tbi" 2>/dev/null || true

# Use -nodownload to force use of the local DB only
java -Xmx4g -jar "$SNPEFF_HOME/snpEff.jar" ann -v -nodownload "$GENOME_ID" "$INPUT_VCF" \
  | bgzip -c > "$OUTPUT_VCF"

# Index if tabix is available
if command -v tabix >/dev/null 2>&1; then
  tabix -p vcf "$OUTPUT_VCF"
fi

echo "Done! Annotated VCF: $OUTPUT_VCF"
