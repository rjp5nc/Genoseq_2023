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
#SBATCH --array=1-404

cd /scratch/rjp5nc/UK2022_2024/whatshap

# Load modules or conda
module load bcftools tabix   # optional if needed
# module load whatshap        # use if available
# Or conda: source activate your_env

# ===== CONFIG =====
CSV="/scratch/rjp5nc/UK2022_2024/daphnia_phylo/Sample_ID_Species_merged_20250627.csv"
VCF="/scratch/rjp5nc/UK2022_2024/daphnia_phylo/trimmed_10bp_repeatmasked_vcf/trimmed10bp_masked_usobtusa.vcf.gz"
REF="/scratch/rjp5nc/Reference_genomes/post_kraken/US_obtusa_onlydaps.fa"
BAM_DIR="/scratch/rjp5nc/UK2022_2024/final_bam_rg2/"
OUT_DIR="/scratch/rjp5nc/UK2022_2024/whatshap/phased_whathap"

mkdir -p $OUT_DIR

# Load the sample for this array task
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" all_samples.txt)

# Find the BAM
BAM=$(ls $BAM_DIR/${SAMPLE}*finalmap_RG.bam 2>/dev/null | head -n1)
OUT="$OUT_DIR/phased_${SAMPLE}.vcf"

if [ -f "$BAM" ]; then
    echo "Phasing $SAMPLE -> $BAM ..."
    whatshap phase \
      --reference $REF \
      --output $OUT \
      $VCF \
      $BAM
    echo "Done phasing $SAMPLE -> $OUT"
else
    echo "Warning: BAM not found for $SAMPLE, skipping..."
fi