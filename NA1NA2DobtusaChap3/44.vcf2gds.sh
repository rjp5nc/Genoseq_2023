#!/usr/bin/env bash
#
#SBATCH -J vcf2gds # A single job name for the array
#SBATCH --ntasks-per-node=10 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-10:00  ### 48 hours
#SBATCH --mem 60G
#SBATCH -o /scratch/rjp5nc/vcf2gds.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/vcf2gds.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

### sbatch ~/Genoseq_2023/vcf2gds.sh
### sacct -j 22867938

module load gcc openmpi R/4.3.1

Rscript --vanilla /home/rjp5nc/Genoseq_2023/finalscripts_manuscript_3/44.vcf2gds.R \
/scratch/rjp5nc/UK2022_2024/buscoanalysis/buscolifted_singlecopy.vcf.gz