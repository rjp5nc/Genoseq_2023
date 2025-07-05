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

#!/bin/bash

module load gatk/4.6.0.0

parameterFile=/scratch/rjp5nc/UK2022_2024/touseforDBI_mito_fullref.csv
wd="/scratch/rjp5nc/UK2022_2024/mitogvcf"

#dos2unix "$parameterFile"
  
#awk -F',' 'NR>1 {print $1 "\t/scratch/rjp5nc/UK2022_2024/mitogvcf/gvcf/" $6 "/" $1 ".g.vcf.gz"}' /scratch/rjp5nc/UK2022_2024/touseforDBI_mito_fullref.csv > /scratch/rjp5nc/UK2022_2024/sample_map.txt


# Extract sample name

#changed these for pulex file

OUTPUT_DIR="/scratch/rjp5nc/UK2022_2024/mito_vcf/merged_gvcf"
GROUP_NAME="dambigua_mito"
REF_FASTA="/scratch/rjp5nc/Reference_genomes/mito_reference/${GROUP_NAME}.fasta"
GVCF_DIR="/scratch/rjp5nc/UK2022_2024/mitogvcf/gvcf/${GROUP_NAME}"
OUT_GVCF="${OUTPUT_DIR}/${GROUP_NAME}_combined.g.vcf.gz"

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Run CombineGVCFs using all .g.vcf.gz files in the folder
gatk --java-options "-Xmx40G" CombineGVCFs \
  -R "$REF_FASTA" \
  $(for f in "$GVCF_DIR"/*.g.vcf.gz; do echo "-V $f"; done) \
  -O "$OUT_GVCF"