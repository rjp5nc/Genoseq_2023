#!/usr/bin/env bash
#
#SBATCH -J bgzip # A single job name for the array
#SBATCH --ntasks-per-node=10 # multi core
#SBATCH -N 1 # on one node
#SBATCH -t 1-00:00 # 1 days
#SBATCH --mem 60G
#SBATCH -o /scratch/rjp5nc/err/m%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/err/m%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

module load bcftools
module load samtools
bgzip /scratch/rjp5nc/UK2022_2024/daphnia_phylo/trimmed_10bp_repeatmasked_vcf/trimmed10bp_masked_usobtusa2.vcf