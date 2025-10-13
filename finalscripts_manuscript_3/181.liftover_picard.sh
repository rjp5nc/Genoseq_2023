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

#species=uspulex
#species=usobtusa
#species=usambigua
species=usobtusa
speciescrossmap=us_obtusa
#vcffile=trimmed10bp_masked_${species}.vcf.gz
#Daphnia_ambigua_Q001_genome.fa
#US_obtusa_onlydaps.fa
#us_pulex_ref_kap4.fa

# Parameters
JAVAMEM=100G
CPU=1

# Liftover using picard
java "-Xmx${JAVAMEM}" -jar $EBROOTPICARD/picard.jar  LiftoverVcf \
I=/scratch/rjp5nc/UK2022_2024/daphnia_phylo/eudobtusa_indv/trimmed10bp_masked_euobtusa.vcf.gz \
O=/scratch/rjp5nc/UK2022_2024/daphnia_phylo/eudobtusa_indv/trimmed10bp_masked_euobtusa_lifted_to_US_picard.vcf.gz \
CHAIN=/scratch/rjp5nc/lastz/eu_obtusa_to_usobtusa/chainnet/liftover.chain \
REJECT=/scratch/rjp5nc/UK2022_2024/daphnia_phylo/eudobtusa_indv/trimmed10bp_masked_euobtusa_lifted_to_US_picard_rejects.vcf.gz \
R=/scratch/rjp5nc/Reference_genomes/post_kraken/assembly.hap2_onlydaps.fasta \
MAX_RECORDS_IN_RAM=10000000 \
WARN_ON_MISSING_CONTIG=true \
RECOVER_SWAPPED_REF_ALT=true 

tabix -p vcf /scratch/rjp5nc/UK2022_2024/daphnia_phylo/trimmed_10bp_repeatmasked_vcf/lifted_${species}.vcf.gz
