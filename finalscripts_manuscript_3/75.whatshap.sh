#!/usr/bin/env bash
#
#SBATCH -J whatshap # A single job name for the array
#SBATCH --ntasks-per-node=10 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-72:00  ### 48 hours
#SBATCH --mem 60G
#SBATCH -o /scratch/rjp5nc/err/whatshap.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/err/whatshap.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

cd /scratch/rjp5nc/UK2022_2024/whatshap

CSV="/scratch/rjp5nc/UK2022_2024/daphnia_phylo/Sample_ID_Species_merged_20250627.csv"
VCF="/scratch/rjp5nc/UK2022_2024/daphnia_phylo/trimmed_10bp_repeatmasked_vcf/trimmed10bp_masked_usobtusa.vcf.gz"
REF="/scratch/rjp5nc/Reference_genomes/post_kraken/US_obtusa_onlydaps.fa"
BAM_DIR="/scratch/rjp5nc/UK2022_2024/final_bam_rg2/"
OUT_DIR="phased_whathap"
MERGED_VCF="$OUT_DIR/phased_merged.vcf.gz"

mkdir -p $OUT_DIR

# Skip header, select SampleID where Species=="Daphnia obtusa" and Continent=="NorthAmerica"
SAMPLES=$(awk -F',' 'NR>1 && $5 ~ /Daphnia obtusa/ && $9 ~ /NorthAmerica/ {gsub(/"/,"",$3); print $3}' $CSV)
SAMPLES2=$(awk -F',' 'NR>1 && $5 ~ /Daphnia obtusa/ && $9 ~ /NorthAmerica/ {gsub(/"/,"",$10); print $10}' $CSV)

ALL_SAMPLES=$(echo -e "$SAMPLES\n$SAMPLES2" | sort | uniq)

echo "Selected samples:"
echo "$ALL_SAMPLES"

# ===== Step 2: Phase each sample =====
for SAMPLE in $ALL_SAMPLES; do
    BAM=$(ls $BAM_DIR/${SAMPLE}*finalmap_RG.bam 2>/dev/null)
    OUT="$OUT_DIR/phased_${SAMPLE}.vcf"
    
    if [ -f "$BAM" ]; then
        echo "Phasing $SAMPLE..."
        whatshap phase \
          --reference $REF \
          --output $OUT \
          $VCF \
          $BAM
        echo "Done phasing $SAMPLE -> $OUT"
    else
        echo "Warning: BAM not found for $SAMPLE, skipping..."
    fi
done

# ===== Step 3: Merge all per-sample phased VCFs =====
echo "Merging all per-sample phased VCFs..."
bcftools merge $OUT_DIR/phased_*.vcf -O z -o $MERGED_VCF
tabix -p vcf $MERGED_VCF

echo "Merged phased VCF ready: $MERGED_VCF"