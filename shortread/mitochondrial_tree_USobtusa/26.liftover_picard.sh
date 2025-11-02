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


#species=dambigua_mito
species=eudobtusa_mito_reverse
#species=kap4Dpulex_mito
species=usdobtusa_mito


# Parameters
JAVAMEM=40G
CPU=1

# Liftover using picard
java "-Xmx${JAVAMEM}" -jar $EBROOTPICARD/picard.jar  LiftoverVcf \
I= /scratch/rjp5nc/UK2022_2024/mito_vcf/merged_vcfs/${species}_genotyped.vcf.gz \
O= /scratch/rjp5nc/UK2022_2024/mito_vcf/lifted_vcfs/lifted_${species}.vcf.gz \
CHAIN=/scratch/rjp5nc/Reference_genomes/mito_reference/lastz/${species}/chainnet/liftover.chain \
REJECT=/scratch/rjp5nc/UK2022_2024/mito_vcf/lifted_vcfs/lifted_${species}_rejected.vcf.gz \
R=/scratch/rjp5nc/Reference_genomes/mito_reference/eudpulex_mito.fasta \
MAX_RECORDS_IN_RAM=10000000 \
WARN_ON_MISSING_CONTIG=true \
RECOVER_SWAPPED_REF_ALT=true 

