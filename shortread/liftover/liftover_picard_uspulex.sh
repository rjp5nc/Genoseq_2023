#!/usr/bin/env bash

#SBATCH -J liftover_picard # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 24:00:00 ### 15 seconds
#SBATCH --mem 110G
#SBATCH -o /scratch/rjp5nc/erroroutputs/nFlo_1.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/erroroutputs/nFlo_1.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

module load tabix
module load picard

#finished eu_obtusa, ambigua, usobtusa

species=uspulex
speciescrossmap=us_pulex
vcffile=trimmed10bp_masked_${species}.vcf.gz
sourcegenome=totalHiCwithallbestgapclosed.fa

# Parameters
JAVAMEM=100G
CPU=1

#!/bin/bash
set -e

contig=$1
input_vcf="/scratch/rjp5nc/UK2022_2024/daphnia_phylo/trimmed_10bp_repeatmasked_vcf/uspulex/trimmed10bp_masked_uspulex.$contig.vcf.gz"
output_vcf="/scratch/rjp5nc/UK2022_2024/daphnia_phylo/trimmed_10bp_repeatmasked_vcf/uspulex/lifted_uspulex.$contig.vcf.gz"
reject_vcf="/scratch/rjp5nc/UK2022_2024/daphnia_phylo/trimmed_10bp_repeatmasked_vcf/uspulex/rejected_uspulex.$contig.vcf.gz"

cd /scratch/rjp5nc/UK2022_2024/daphnia_phylo/trimmed_10bp_repeatmasked_vcf/uspulex/

java "-Xmx${JAVAMEM}" -jar $EBROOTPICARD/picard.jar LiftoverVcf \
  I=$input_vcf \
  O=$output_vcf \
  CHAIN=/scratch/rjp5nc/lastz/$speciescrossmap/chainnet/liftover.chain \
  REJECT=$reject_vcf \
  R=/scratch/rjp5nc/Reference_genomes/post_kraken/$sourcegenome \
  MAX_RECORDS_IN_RAM=10000000 \
  WARN_ON_MISSING_CONTIG=true \
  RECOVER_SWAPPED_REF_ALT=true

tabix -p vcf $output_vcf