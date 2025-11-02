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

awk 'NR>3 {print $5"\t"$6-1"\t"$7}' /scratch/rjp5nc/removedups/eu_dobtusa/assembly.hap2_onlydaps.fasta.out > /scratch/rjp5nc/removedups/eu_dobtusa/euobtusa_repeats.bed

bedtools subtract -header \
  -a /scratch/rjp5nc/UK2022_2024/daphnia_phylo/eudobtusa_indv/trimmed10bp_allsites_euobtusa.vcf.gz \
  -b /scratch/rjp5nc/removedups/eu_dobtusa/euobtusa_repeats.bed \
| bgzip > /scratch/rjp5nc/UK2022_2024/daphnia_phylo/eudobtusa_indv/trimmed10bp_masked_euobtusa.vcf.gz

tabix -p vcf /scratch/rjp5nc/UK2022_2024/daphnia_phylo/eudobtusa_indv/trimmed10bp_masked_euobtusa.vcf.gz
