#!/usr/bin/env bash

#SBATCH -J windowed_het # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 3-0:00:00 ### 15 seconds
#SBATCH --mem 60G
#SBATCH -o /scratch/rjp5nc/erroroutputs/beast.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/erroroutputs/beast.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab
#SBATCH --mail-type=END               # Send email at job completion
#SBATCH --mail-user=rjp5nc@virginia.edu    # Email address for notifications

module load gatk
cd /scratch/rjp5nc/UK2022_2024/daphnia_phylo/eupulex_indv/gvcf


find /scratch/rjp5nc/UK2022_2024/daphnia_phylo/gvcf/eupulex_chr/Scaffold_* \
    -name "*.g.vcf.gz" > all_gvcfs.list

gatk CombineGVCFs \
    -R /scratch/rjp5nc/Reference_genomes/orig_ref/eu_pulex_totalHiCwithallbestgapclosed.clean.fa \
    $(sed 's/^/-V /' all_gvcfs.list) \
    -O eupulex_all_scaffolds.merged.g.vcf.gz
