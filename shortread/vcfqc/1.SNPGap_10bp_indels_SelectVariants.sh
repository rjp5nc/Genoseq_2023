#!/usr/bin/env bash
#
#SBATCH -J SNP_filtering
#SBATCH --ntasks-per-node=4 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-72:00 # hours
#SBATCH --mem 10G
#SBATCH -o /scratch/rjp5nc/err/FilterVCFs.out # Standard output
#SBATCH -e /scratch/rjp5nc/err/FilterVCFs.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab
#SBATCH --array=1-369%50
#SBATCH --mail-type=END               # Send email at job completion
#SBATCH --mail-user=rjp5nc@virginia.edu    # Email address for notifications


# This script will remove SNPS within 10 base pairs of indels by a chromosome basis using bcftools
# Also filters out only SNPs using gatk

# Modules to load
module load bcftools/1.17
module load gatk/4.6.0.0

# Working directory

species=usobtusa_vcf
ref=/scratch/rjp5nc/Reference_genomes/post_kraken/US_obtusa_onlydaps.fa

wd=/scratch/rjp5nc/UK2022_2024/daphnia_phylo/10bp_vcf/$species

mkdir -p /scratch/rjp5nc/UK2022_2024/daphnia_phylo/10bp_vcf/$species

WORKING_FOLDER=/scratch/rjp5nc/UK2022_2024/daphnia_phylo/vcf/$species
# Intervals to analyze
intervals="/scratch/rjp5nc/UK2022_2024/daphnia_phylo/interval_DBI_paramList_usobtusa.txt"

# Parameters
JAVAMEM=10G
CPU=4


# Chromosome
chrom=$( cat ${intervals} | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f2 )

# Start
start=$( cat ${intervals} | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f3 )

# Stop
stop=$( cat ${intervals} | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f4 )

# Input file
IN_GZVCF=$WORKING_FOLDER/${chrom}.${start}.${stop}.vcf.gz

# Move to working directory

cd ${wd}

# Remove SNPs 10 bp from indels
bcftools filter \
--SnpGap 10 \
${IN_GZVCF} \
--output ${wd}/${chrom}.${start}.${stop}_filtsnps10bpindels.vcf.gz \
-O z \
--threads ${CPU}

module load gcc/14.2.0 htslib/1.17

# Index filtered vcf
tabix -p vcf ${wd}/${chrom}.${start}.${stop}_filtsnps10bpindels.vcf.gz

mkdir -p ${wd}/${species}_snps

# Filter only SNPs
gatk --java-options "-Xmx${JAVAMEM}" SelectVariants \
-V ${wd}/${chrom}.${start}.${stop}_filtsnps10bpindels.vcf.gz \
--select-type-to-include SNP \
-R $ref \
-O ${wd}/${species}_snps/${chrom}.${start}.${stop}_filtsnps10bpindels_snps.vcf.gz