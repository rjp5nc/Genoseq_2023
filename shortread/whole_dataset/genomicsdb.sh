#!/usr/bin/env bash
#
#SBATCH -J genomicsdb # A single job name for the array
#SBATCH --ntasks-per-node=10 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 5-12:00:00 # 8 hours
#SBATCH --mem 80G
#SBATCH -o /scratch/rjp5nc/err/gatk.chrom.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/err/gatk.chrom.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab
#SBATCH --array=1-113%50


# This script will merge gVCFs into a unified database for genotype calling.
# This will be done using a per chromosome approach

# Load modules
module load gatk/4.6.0.0

# Name of pipeline
PIPELINE="GenomicsDBImport"

# Working folder is core folder where this pipeline is being run.
WORKING_FOLDER="/scratch/rjp5nc/UK2022_2024/daphnia_phylo/DBI_usobtusa"

#cat -A /scratch/rjp5nc/UK2022_2024/daphnia_phylo/gvcf/euobtusa_chr/euobtusa_gvcflist.txt
#sed 's/\r$//' /scratch/rjp5nc/UK2022_2024/daphnia_phylo/gvcf/euobtusa_chr/euobtusa_gvcflist.txt | sed 's/\s*$//' > /scratch/rjp5nc/UK2022_2024/daphnia_phylo/gvcf/euobtusa_chr/euobtusa_gvcflist.clean.txt
#cat -A /scratch/rjp5nc/UK2022_2024/daphnia_phylo/gvcf/euobtusa_chr/euobtusa_gvcflist.clean.txt | head -n 20
#tr '\t' '   ' < /scratch/rjp5nc/UK2022_2024/daphnia_phylo/gvcf/euobtusa_chr/euobtusa_gvcflist.clean.txt > /scratch/rjp5nc/UK2022_2024/daphnia_phylo/gvcf/euobtusa_chr/euobtusa_gvcflist.clean2.txt
#sed 's/ *$//' /scratch/rjp5nc/UK2022_2024/daphnia_phylo/gvcf/euobtusa_chr/euobtusa_gvcflist.txt > /scratch/rjp5nc/UK2022_2024/daphnia_phylo/gvcf/euobtusa_chr/euobtusa_gvcflist2.txt

# Chromosomes to analyze
intervals="/scratch/rjp5nc/UK2022_2024/daphnia_phylo/interval_DBI_paramList_usobtusa.txt"

# Parameters
JAVAMEM=80G
CPU=10

# Move to working directory
cd $WORKING_FOLDER

# Chromosome
i=$( cat ${intervals} | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f2 )

# Start
start=$( cat ${intervals} | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f3 )

# Stop
stop=$( cat ${intervals} | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f4 )

if [[ -d "TEMP_Daphnia_DBI_${i}_${start}_${stop}" ]]
then
echo "Working TEMP_Daphnia_DBI folder exist"
echo "lets move on"
date
else
echo "folder doesnt exist. lets fix that"
mkdir $WORKING_FOLDER/TEMP_Daphnia_DBI_${i}_${start}_${stop}
date
fi

echo ${i}:${start}-${stop} "is being processed" $(date)

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
--genomicsdb-workspace-path $WORKING_FOLDER/Daphnia_DBI_${i}_${start}_${stop} \
--tmp-dir $WORKING_FOLDER/TEMP_Daphnia_DBI_${i}_${start}_${stop} \
--batch-size 50 \
--sample-name-map /scratch/rjp5nc/UK2022_2024/daphnia_phylo/samplemapnames/${i}.txt \
--reader-threads $CPU \
-L ${i}:${start}-${stop}

#sample name from g.vcf.gz 

# Remove temp workspace
rm -rf $WORKING_FOLDER/TEMP_Daphnia_DBI_${i}_${start}_${stop}

echo ${i} "done" $(date)