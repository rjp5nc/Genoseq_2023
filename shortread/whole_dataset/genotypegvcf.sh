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
#SBATCH --array=1-254%50
#SBATCH --mail-type=END               # Send email at job completion
#SBATCH --mail-user=rjp5nc@virginia.edu    # Email address for notifications

# This script will conduct genotype calling on the GenomeDBI object

# Load Modules
module load gatk/4.6.0.0

# Parameters
JAVAMEM=80G
CPU=10

# Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/scratch/rjp5nc/UK2022_2024/daphnia_phylo/vcf

# Reference genome
REFERENCE=/scratch/rjp5nc/Reference_genomes/post_kraken/assembly.hap2_onlydaps.fasta

# Intervals to analyze
intervals="/scratch/rjp5nc/UK2022_2024/daphnia_phylo/interval_DBI_paramList_euobtusa.txt"

species="DBI_euobtusa"
# This part of the pipeline will generate log files to record warnings and completion status

# Move to working directory
cd $WORKING_FOLDER

# Chromosome
i=$( cat ${intervals} | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f2 )

# Start
start=$( cat ${intervals} | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f3 )

# Stop
stop=$( cat ${intervals} | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f4 )

# Create temp folder
if [[ -d "TEMP_Daphnia_Genotype_${i}_${start}_${stop}" ]]
then
echo "Working TEMP_Daphnia_Genotype folder exist"
echo "Lets move on"
date
else
echo "Folder doesnt exist. Lets fix that"
mkdir $WORKING_FOLDER/TEMP_Daphnia_Genotype_${i}_${start}_${stop}
date
fi

echo ${i}_${start}_${stop} "is being processed" $(date)

# Identify the Genome database to genotyoe
GenomeDB_path=`echo /scratch/rjp5nc/UK2022_2024/daphnia_phylo/$species/Daphnia_DBI_${i}_${start}_${stop}`

# Genotype call the samples in the DBI merged set
gatk --java-options "-Xmx${JAVAMEM}" GenotypeGVCFs \
-R $REFERENCE \
-V gendb://$GenomeDB_path \
--tmp-dir $WORKING_FOLDER/TEMP_Daphnia_Genotype_${i}_${start}_${stop} \
-O $WORKING_FOLDER/${i}.${start}.${stop}.vcf.gz \
#--genomicsdb-use-vcf-codec \
-L ${i}:${start}-${stop}

# Remove temp folder
rm -rf $WORKING_FOLDER/TEMP_Daphnia_Genotype_${i}_${start}_${stop}

echo ${i} "done" $(date)


#find . -type f -name "*.vcf.gz" | sed 's|^\./||' > /scratch/rjp5nc/UK2022_2024/daphnia_phylo/unmerged_eudobtusa_vcf_files.txt

