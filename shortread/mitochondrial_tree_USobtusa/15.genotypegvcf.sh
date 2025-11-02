#!/usr/bin/env bash
#
#SBATCH -J genotypegvcf # A single job name for the array
#SBATCH --ntasks-per-node=10 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 7:00:00 # 8 hours
#SBATCH --mem 80G
#SBATCH -o /scratch/rjp5nc/err/genotypegvcf.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/err/genotypegvcf.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab
#SBATCH --mail-type=END               # Send email at job completion
#SBATCH --mail-user=rjp5nc@virginia.edu    # Email address for notifications

# This script will conduct genotype calling on the GenomeDBI object

# Load Modules
module load gatk/4.6.0.0

# Parameters
JAVAMEM=40G
CPU=10

#NEED TO DO USPULEX/AMBIGUA

# Working folder is core folder where this pipeline is being run.
folder=usdobtusa_mito

  gatk --java-options "-Xmx40G" GenotypeGVCFs \
    -R "/scratch/rjp5nc/Reference_genomes/mito_reference/${folder}.fasta" \
    -V "/scratch/rjp5nc/UK2022_2024/mito_vcf/merged_gvcf/${folder}_combined.g.vcf.gz" \
    -O "/scratch/rjp5nc/UK2022_2024/mito_vcf/merged_vcfs/${folder}_genotyped.vcf.gz"


#Flag to spit out all sites?