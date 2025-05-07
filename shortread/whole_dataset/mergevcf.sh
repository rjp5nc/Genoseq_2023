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
WORKING_FOLDER=/scratch/rjp5nc/UK2022_2024/daphnia_phylo/vcf

# Intervals to analyze
# Create VCF combine list
ls $WORKING_FOLDER/*vcf.gz > $WORKING_FOLDER/interval_paramList.list

# VCF list location
intervals=/scratch/rjp5nc/UK2022_2024/daphnia_phylo/unmerged_eudobtusa_vcf_files.txt

# Combined VCF name
output="/scratch/rjp5nc/UK2022_2024/daphnia_phylo/raw_vcf/raw_euobtusa.vcf.gz"

# Gatk parameters
JAVAMEM=100G
CPU=1

# Move to working directory
cd $WORKING_FOLDER

# Combine VCFs using MergeVcfs
gatk --java-options "-Xmx${JAVAMEM}" MergeVcfs \
-I $intervals \
-O $output

# Finish
echo "Complete" $(date)
