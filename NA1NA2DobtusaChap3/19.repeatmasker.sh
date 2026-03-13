#!/usr/bin/env bash
#
#SBATCH -J RemoveRepeats # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-2:00  ### 1 hours
#SBATCH --mem 24G
#SBATCH -o /scratch/rjp5nc/UK2022_2024/allshortreads/chr/2022seq.concat.Removereps.renamed.vcf.gz
#SBATCH -e /scratch/rjp5nc/outputerrors/download_%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

### sbatch ~/Genoseq_2023/RemoveRepetitiveRegions.sh
### sacct -j 22867938
### scancel
### cat

module load gcc/12.4.0 bedtools/2.30.0 htslib

awk 'NR>3 {print $5"\t"$6-1"\t"$7}' /scratch/rjp5nc/removedups/us_dobtusa/US_obtusa_onlydaps.fa.out > /scratch/rjp5nc/removedups/us_dobtusa/usobtusa_repeats.bed

bedtools subtract -header \
-a /project/berglandlab/Robert/UKSequencing2022_2024/trimmed10bp_onlysnp_vcf/trimmed10bp_usobtusa_vcf.vcf.gz \
-b /scratch/rjp5nc/removedups/us_dobtusa/usobtusa_repeats.bed > /scratch/rjp5nc/UK2022_2024/daphnia_phylo/trimmed_10bp_repeatmasked_vcf/trimmed10bp_masked_usobtusa.vcf.gz


tabix -p vcf  /scratch/rjp5nc/UK2022_2024/daphnia_phylo/trimmed_10bp_repeatmasked_vcf/trimmed10bp_masked_usobtusa.vcf.gz


bgzip -c /scratch/rjp5nc/UK2022_2024/daphnia_phylo/trimmed_10bp_repeatmasked_vcf/trimmed10bp_masked_usobtusa.vcf.gz \
  > /scratch/rjp5nc/UK2022_2024/daphnia_phylo/trimmed_10bp_repeatmasked_vcf/trimmed10bp_masked_usobtusa.bgz.vcf.gz

# Index with tabix
tabix -p vcf /scratch/rjp5nc/UK2022_2024/daphnia_phylo/trimmed_10bp_repeatmasked_vcf/trimmed10bp_masked_usobtusa.bgz.vcf.gz