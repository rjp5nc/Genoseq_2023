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

# Parameters

# Move to working directory
cd $WORKING_FOLDER

# Chromosome

WORKING_FOLDER="/scratch/rjp5nc/UK2022_2024/daphnia_phylo/DBI_mitos/${folder}"

mkdir -p "${WORKING_FOLDER}" 
mkdir -p "$WORKING_FOLDER/temp"
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


sample_list="/scratch/rjp5nc/UK2022_2024/sample_map_fixed.txt"

# Read line and split into two variables
read -r sample_id gvcf_path < <(awk -F'\t' -v line=${SLURM_ARRAY_TASK_ID} 'NR==line' ${sample_list})

echo "Sample ID: $sample_id"
echo "GVCF Path: $gvcf_path"

interval=$(awk '{printf "%s:%d-%d", $1, $2+1, $3}' /scratch/rjp5nc/Reference_genomes/mito_reference/${folder}.bed)
echo "$interval"

mkdir -p "${WORKING_FOLDER}/TEMP_Daphnia_DBI_${species}_${sample_id}"

gatk --java-options "-Xmx${JAVAMEM}" GenomicsDBImport \
--genomicsdb-workspace-path $WORKING_FOLDER/${species}/Daphnia_DBI_${species} \
--tmp-dir $WORKING_FOLDER/TEMP_Daphnia_DBI_${species}_${sample_id} \
--batch-size 50 \
--sample-name-map /scratch/rjp5nc/UK2022_2024/sample_map_fixed.txt \
--reader-threads $CPU \
-L "${interval}"

#sample name from g.vcf.gz 
# Remove temp workspace
rm -rf $WORKING_FOLDER/TEMP_Daphnia_DBI_${i}_${start}_${stop}

echo ${i} "done" $(date)