#!/usr/bin/env bash

#SBATCH -J liftover_picard # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 1:00:00 ### 15 seconds
#SBATCH --mem 100G
#SBATCH -o /scratch/rjp5nc/erroroutputs/nFlo_1.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/erroroutputs/nFlo_1.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

module load tabix
module load picard


#species=uspulex
#species=usobtusa
#species=usambigua
species=euobtusa
speciescrossmap=eu_obtusa
vcffile=trimmed10bp_masked_${species}.vcf.gz
sourcegenome=totalHiCwithallbestgapclosed.fa
#Daphnia_ambigua_Q001_genome.fa
#US_obtusa_onlydaps.fa
#us_pulex_ref_kap4.fa

# Parameters
JAVAMEM=100G
CPU=1

# Liftover using picard
java "-Xmx${JAVAMEM}" -jar $EBROOTPICARD/picard.jar  LiftoverVcf \
I=/scratch/rjp5nc/UK2022_2024/daphnia_phylo/trimmed_10bp_repeatmasked_vcf/$vcffile \
O=/scratch/rjp5nc/UK2022_2024/daphnia_phylo/trimmed_10bp_repeatmasked_vcf/lifted_${species}.vcf.gz \
CHAIN=/scratch/rjp5nc/lastz/$speciescrossmap/chainnet/liftover.chain \
REJECT=/scratch/rjp5nc/UK2022_2024/daphnia_phylo/trimmed_10bp_repeatmasked_vcf/lifted_${species}_rejected.vcf \
R=/scratch/rjp5nc/Reference_genomes/post_kraken/$sourcegenome \
MAX_RECORDS_IN_RAM=50000000 \
WARN_ON_MISSING_CONTIG=true \
RECOVER_SWAPPED_REF_ALT=true 

tabix -p vcf /scratch/rjp5nc/UK2022_2024/daphnia_phylo/trimmed_10bp_repeatmasked_vcf/lifted_${species}.vcf.gz
