#!/usr/bin/env bash

#SBATCH -J merge2 # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 24:00:00 ### 15 seconds
#SBATCH --mem 50G
#SBATCH -o /scratch/rjp5nc/erroroutputs/nFlo_1.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/erroroutputs/nFlo_1.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

module load gatk/4.6.0.0  # adjust version as needed

species=uspulex
speciescrossmap=us_pulex
vcffile=trimmed10bp_masked_${species}.vcf.gz
sourcegenome=totalHiCwithallbestgapclosed.fa

cd /scratch/rjp5nc/UK2022_2024/daphnia_phylo/trimmed_10bp_repeatmasked_vcf/lifted_vcf

# Define species name
input_dir="."  # or wherever your VCFs are
output_vcf="lifted_all.vcf.gz"

# Create a list of input files
inputs=$(ls ${input_dir}/*.vcf.gz | sed 's/^/--INPUT /' | tr '\n' ' ')

# Merge
gatk MergeVcfs \
  ${inputs} \
  --OUTPUT ${output_vcf}

# Index
gatk IndexFeatureFile -I ${output_vcf}