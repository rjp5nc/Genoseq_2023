#!/usr/bin/env bash

#SBATCH -J BEAST # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 6-0:00:00 ### 15 seconds
#SBATCH --mem 60G
#SBATCH -o /scratch/rjp5nc/erroroutputs/beast.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/erroroutputs/beast.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab
#SBATCH --mail-type=END               # Send email at job completion
#SBATCH --mail-user=rjp5nc@virginia.edu    # Email address for notifications



# Optional: create a new environment
conda create -n snpeff_env openjdk=11 -y

# Activate it
conda activate snpeff_env
conda install -c bioconda snpeff -y
snpEff -h



mkdir -p /scratch/rjp5nc/snpEff/data/US_obtusa

cp /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/Daphnia_obtusa_genes/Daphnia_obtusa_FS6_genome.fa /scratch/rjp5nc/snpEff/data/US_obtusa/sequences.fa
cp /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/Daphnia_obtusa_genes/Daphnia_obtusa_FS6_genome.gtf /scratch/rjp5nc/snpEff/data/US_obtusa/genes.gtf




~/miniconda3/envs/snpeff_env/bin/snpEff \
-c /scratch/rjp5nc/snpEff/snpEff.config US_obtusa /path/to/input.vcf > /scratch/rjp5nc/snpEff/output.ann.vcf
