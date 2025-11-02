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

#vcffile=trimmed10bp_euobtusa_vcf.vcf.gz_busco.vcf.gz
vcffile=trimmed10bp_usambigua_vcf.vcf.gz_busco.vcf.gz
#vcffile=trimmed10bp_usobtusa_vcf.vcf.gz_busco.vcf.gz
#vcffile=trimmed10bp_uspulex_vcf.vcf.gz_busco.vcf.gz


#speciescrossmap=eu_obtusa
speciescrossmap=us_ambigua
#speciescrossmap=us_obtusa
#speciescrossmap=us_pulex

sourcegenome=totalHiCwithallbestgapclosed.fa

# Parameters
JAVAMEM=40G
CPU=1

# Liftover using picard
java "-Xmx${JAVAMEM}" -jar $EBROOTPICARD/picard.jar  LiftoverVcf \
I= /scratch/rjp5nc/UK2022_2024/buscoanalysis/filteredvcf/$vcffile \
O= /scratch/rjp5nc/UK2022_2024/buscoanalysis/filteredvcf_lifted/lifted_${speciescrossmap}.vcf.gz \
CHAIN=/scratch/rjp5nc/lastz/$speciescrossmap/chainnet/liftover.chain \
REJECT=/scratch/rjp5nc/UK2022_2024/buscoanalysis/filteredvcf_lifted/lifted_${speciescrossmap}_rejected.vcf.gz \
R=/scratch/rjp5nc/Reference_genomes/post_kraken/$sourcegenome \
MAX_RECORDS_IN_RAM=10000000 \
WARN_ON_MISSING_CONTIG=true \
RECOVER_SWAPPED_REF_ALT=true 

tabix -p vcf /scratch/rjp5nc/UK2022_2024/buscoanalysis/filteredvcf_lifted/lifted_${speciescrossmap}.vcf.gz
