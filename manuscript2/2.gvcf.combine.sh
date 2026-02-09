#!/usr/bin/env bash

#SBATCH -J windowed_het # A single job name for the array
#SBATCH --ntasks-per-node=10 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 3-0:00:00 ### 15 seconds
#SBATCH --mem 100G
#SBATCH -o /scratch/rjp5nc/erroroutputs/beast.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/erroroutputs/beast.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab
#SBATCH --mail-type=END               # Send email at job completion
#SBATCH --mail-user=rjp5nc@virginia.edu    # Email address for notifications
#SBATCH --array=1-12


module load gatk
cd /scratch/rjp5nc/UK2022_2024/daphnia_phylo/eupulex_indv/gvcf

find combined_by_scaffold -name "Scaffold_*.merged.g.vcf.gz" | sort -V > merged_gvcfs.list

gatk GatherVcfs \
  -I merged_gvcfs.list \
  -O eupulex_all_scaffolds.merged.g.vcf.gz

tabix -p vcf eupulex_all_scaffolds.merged.g.vcf.gz
