#!/usr/bin/env bash
#
#SBATCH -J pixy # A single job name for the array
#SBATCH --ntasks-per-node=20 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-92:00  ### 48 hours
#SBATCH --mem 80G
#SBATCH -o /scratch/rjp5nc/err/pixy.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/err/pixy.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab


module load bcftools
cd /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv

gunzip -c trimmed10bp_allsites_usobtusa.vcf.gz | bgzip -@ 16 -c > trimmed10bp_allsites_usobtusa.bgz.vcf.gz
