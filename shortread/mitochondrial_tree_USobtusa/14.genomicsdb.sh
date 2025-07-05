#!/usr/bin/env bash
#
#SBATCH -J genomicsdb # A single job name for the array
#SBATCH --ntasks-per-node=10 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-12:00:00 # 8 hours
#SBATCH --mem 40G
#SBATCH -o /scratch/rjp5nc/err/gatk.chrom.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/err/gatk.chrom.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab
#SBATCH --array=1-521%50


# This script will merge gVCFs into a unified database for genotype calling.
# This will be done using a per chromosome approach


#NEED TO DO AMBIGUA

# Load modules
module load gatk/4.6.0.0

# Name of pipeline
PIPELINE="GenomicsDBImport"

# Working folder is core folder where this pipeline is being run.

#folder=dambigua_mito
#folder=eudobtusa_mito
#folder=eudpulex_mito
#folder=kap4Dpulex_mito
folder=usdobtusa_mito

WORKING_FOLDER="/scratch/rjp5nc/UK2022_2024/daphnia_phylo/DBI_mitos/${folder}"

mkdir -p "${WORKING_FOLDER}" 

# Parameters
JAVAMEM=40G
CPU=10

# Move to working directory
cd $WORKING_FOLDER

#cd /scratch/rjp5nc/UK2022_2024/daphnia_phylo/samplemapnames

#while IFS= read -r line; do
#  # Extract the basename of the directory (e.g., h2tg000008l)
#  chr=$(echo "$line" | awk -F'/' '{for(i=1;i<=NF;i++) if($i ~ /^h2tg[0-9]+[lc]$/) print $i}')

#  # Write the line to the corresponding .txt file
#  echo "$line" >> "${chr}.txt"
#done < "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/gvcf/euobtusapathsfixed.txt"

# Merge VCFs using GenomicsDBImport

#zcat /scratch/rjp5nc/UK2022_2024/daphnia_phylo/gvcf/euobtusa_chr/h2tg000002l/Gilmer5_H9.h2tg000002l.1018082.g.vcf.gz | head -n 500


gatk --java-options "-Xmx${JAVAMEM}" GenomicsDBImport \
  --genomicsdb-workspace-path $WORKING_FOLDER \
  -L /scratch/rjp5nc/Reference_genomes/mito_reference/${folder}.bed \
  --sample-name-map /scratch/rjp5nc/UK2022_2024/sample_map.txt \
  --reader-threads $CPU

#sample name from g.vcf.gz 
# Remove temp workspace
rm -rf $WORKING_FOLDER/TEMP_Daphnia_DBI_${i}_${start}_${stop}

echo ${i} "done" $(date)