#!/usr/bin/env bash
#
#SBATCH -J trim # A single job name for the array
#SBATCH --ntasks-per-node=10 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-10:00  ### 48 hours
#SBATCH --mem 50G
#SBATCH -o /scratch/rjp5nc/err/trim.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/err/trim.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab
#SBATCH --array=1-<NUM_SAMPLES>

#NUM_SAMPLES=$(wc -l < /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/final_vcf_filter_two_of_each.txt)
# submit SLURM array job
#sbatch --array=1-${NUM_SAMPLES} 89.VCF_to_fasta.sh

module load bcftools

#SLURM_ARRAY_TASK_ID=1

VCF="/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/trimmed10bp_filtered_two_of_each.vcf.gz"
REF="/scratch/rjp5nc/Reference_genomes/post_kraken/First_12_US_obtusa_onlydaps.fasta"
SAMPLE_LIST="/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/final_vcf_filter_two_of_each.txt"
OUTDIR="/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/consensus_fastas"

# make sure output directory exists
mkdir -p "$OUTDIR"

# get the sample name corresponding to this array task
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST")

# run bcftools consensus
bcftools consensus -f "$REF" -s "$SAMPLE" -H I "$VCF" > "$OUTDIR/${SAMPLE}.fasta"