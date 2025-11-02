#!/usr/bin/env bash
#
#SBATCH -J popgenome # A single job name for the array
#SBATCH --ntasks-per-node=20 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-100:00  ### 48 hours
#SBATCH --mem 100G
#SBATCH -o /scratch/rjp5nc/err/popgenome.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/err/popgenome.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

### sbatch ~/Genoseq_2023/vcf2gds.sh
### sacct -j 22867938

module load gcc/12.4.0 bedtools/2.30.0 htslib bcftools

cd /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/

bcftools query -l  trimmed10bp_allsites_usobtusa.bgz.vcf.gz > vcf_samples_old.txt



bcftools view -S passed_superclones.txt -Oz \
    --threads 16 \
    -o trimmed10bp_allsites_usobtusa.filtered_bgz.vcf.gz \
    trimmed10bp_allsites_usobtusa.bgz.vcf.gz

#bedtools subtract -header \
#-a /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/trimmed10bp_allsites_usobtusa.filtered_bgz.vcf.gz \
#-b /scratch/rjp5nc/removedups/us_dobtusa/usobtusa_repeats.bed > /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/trimmed10bp_allsites_Repeatmasked_usobtusa.filtered_bgz.vcf.gz

#tabix -p vcf /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/trimmed10bp_allsites_Repeatmasked_usobtusa.filtered_bgz.vcf.gz

#bgzip -c /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/trimmed10bp_allsites_Repeatmasked_usobtusa.filtered_bgz.vcf.gz \
#  > /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/trimmed10bp_allsites_Repeatmasked_usobtusa.filtered_bgz2.vcf.gz



bedtools subtract -header \
  -a /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/trimmed10bp_allsites_usobtusa.filtered_bgz.vcf.gz \
  -b /scratch/rjp5nc/removedups/us_dobtusa/usobtusa_repeats.bed \
| bgzip -c > /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/trimmed10bp_allsites_Repeatmasked_usobtusa.filtered_bgz.vcf.gz

tabix -p vcf /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/trimmed10bp_allsites_Repeatmasked_usobtusa.filtered_bgz.vcf.gz
