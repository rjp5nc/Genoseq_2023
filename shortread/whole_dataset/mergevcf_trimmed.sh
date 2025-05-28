#!/usr/bin/env bash
#
#SBATCH -J CombineVCFs # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 5-00:00
#SBATCH --mem 100G
#SBATCH -o /scratch/rjp5nc/err/genotypegvcf.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/err/genotypegvcf.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab
#SBATCH --mail-type=END               # Send email at job completion
#SBATCH --mail-user=rjp5nc@virginia.edu    # Email address for notifications


# This script will merge all final VCFs together.

# Load modules
module load gatk/4.6.0.0

# Working folder is core folder where this pipeline is being run.

#euobtusa_vcf  eupulex_vcf  eupulex_vcf_old  usambigua_vcf  usobtusa_vcf  uspulex_vcf
#euobtusa_vcf  eupulex_vcf  eupulex_vcf_old  usambigua_vcf  usobtusa_vcf  uspulex_vcf
species=eupulex_vcf_old

WORKING_FOLDER=/scratch/rjp5nc/UK2022_2024/daphnia_phylo/10bp_vcf/$species/${species}_snps


# Intervals to analyze
# Create VCF combine list
cd $WORKING_FOLDER
ls *vcf.gz > $WORKING_FOLDER/${species}interval_paramList.list

# VCF list location
intervals=$WORKING_FOLDER/${species}interval_paramList.list

# Combined VCF name
output="/scratch/rjp5nc/UK2022_2024/daphnia_phylo/trimmed10bp_vcf/trimmed10bp_${species}.vcf.gz"

# Gatk parameters
JAVAMEM=40G
CPU=20

# Move to working directory
cd $WORKING_FOLDER

# Combine VCFs using MergeVcfs
gatk --java-options "-Xmx${JAVAMEM}" MergeVcfs \
-I $intervals \
-O $output

# Finish
echo "Complete" $(date)


