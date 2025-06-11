#!/usr/bin/env bash
#
#SBATCH -J split_pulex
#SBATCH --ntasks-per-node=4 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-72:00 # hours
#SBATCH --mem 40G
#SBATCH -o /scratch/rjp5nc/err/FilterVCFs.out # Standard output
#SBATCH -e /scratch/rjp5nc/err/FilterVCFs.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab
#SBATCH --array=1-13

module load bcftools

cd /scratch/rjp5nc/UK2022_2024/daphnia_phylo/trimmed_10bp_repeatmasked_vcf

VCF=trimmed10bp_masked_uspulex.vcf.gz
OUTDIR=uspulex
CONTIG=$(sed -n "${SLURM_ARRAY_TASK_ID}p" contigs.txt)

mkdir -p "$OUTDIR"

bcftools view -r "$CONTIG" -Oz -o ${OUTDIR}/trimmed10bp_masked_uspulex.${CONTIG}.vcf.gz "$VCF"
