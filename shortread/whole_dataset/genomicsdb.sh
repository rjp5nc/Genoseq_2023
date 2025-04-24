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

# This script will merge gVCFs into a unified database for genotype calling.
# This will be done using a per chromosome approach





find /scratch/rjp5nc/UK2022_2024/daphnia_phylo/gvcf/euobtusa_chr -type f -name "*.g.vcf.gz" | while read file; do
    printf "%s\t%s\n" "$(basename "$file")" "$file"
done > 


find /scratch/rjp5nc/UK2022_2024/daphnia_phylo/gvcf -type f -name "*.g.vcf.gz" | \
awk '{ file=$0; split(file, a, "/"); printf "%-40s %s\n", a[length(a)], file }' > gvcf_euobtusa_list.txt

find /scratch/rjp5nc/UK2022_2024/daphnia_phylo/gvcf -type f -name "*.g.vcf.gz" | \
while IFS= read -r filepath; do
    filename=$(basename "$filepath")
    echo -e "$filename\t$filepath"
done > gvcf_file_list.txt


# Load modules
module load gatk/4.1.6.0

# Name of pipeline
PIPELINE="GenomicsDBImport"

# Working folder is core folder where this pipeline is being run.
WORKING_FOLDER="/scratch/rjp5nc/UK2022_2024/daphnia_phylo/DBI"

# Chromosomes to analyze
intervals=/scratch/rjp5nc/UK2022_2024/daphnia_phylo/interval_DBI_paramList

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

# Merge VCFs using GenomicsDBImport
gatk --java-options "-Xmx${JAVAMEM}" GenomicsDBImport \
--genomicsdb-workspace-path $WORKING_FOLDER/Daphnia_DBI_${i}_${start}_${stop} \
--tmp-dir $WORKING_FOLDER/TEMP_Daphnia_DBI_${i}_${start}_${stop} \
--batch-size 50 \
--sample-name-map /scratch/csm6hg/daphnia_phylo/data/interval.total.test.txt \
--reader-threads $CPU \
-L ${i}:${start}-${stop}

# Remove temp workspace
rm -rf $WORKING_FOLDER/TEMP_Daphnia_DBI_${i}_${start}_${stop}

echo ${i} "done" $(date)