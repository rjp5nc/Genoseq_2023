#!/usr/bin/env bash
#
#SBATCH -J Busco # A single job name for the array
#SBATCH --ntasks-per-node=10 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-05:00 # 5 hours
#SBATCH --mem 60G
#SBATCH -o /scratch/rjp5nc/err/busco.out # Standard output
#SBATCH -e /scratch/rjp5nc/err/busco.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

#conda create -n busco_env -c bioconda -c conda-forge busco
conda activate busco_env
# apptainer pull docker://ezlabgva/busco:v5.4.7_cv1
module load apptainer
module load bcftools
# Move to data directory

cd /scratch/rjp5nc/UK2022_2024/buscoanalysis

#Ref=US_obtusa_onlydaps.fa
#Ref=totalHiCwithallbestgapclosed.fa
Ref=Daphnia_ambigua_Q001_genome.fa
#Ref=assembly.hap2_onlydaps.fasta
#Ref=us_pulex_ref_kap4.fa

#vcf=trimmed10bp_usobtusa_vcf.vcf.gz
#vcf=trimmed10bp_eupulex_vcf_old.vcf.gz
#vcf=trimmed10bp_eupulex_vcf.vcf.gz 
vcf=trimmed10bp_usambigua_vcf.vcf.gz
#vcf=trimmed10bp_euobtusa_vcf.vcf.gz
#vcf=trimmed10bp_uspulex_vcf.vcf.gz

bcftools view -R /scratch/rjp5nc/UK2022_2024/buscoanalysis/${Ref}_BUSCO_combined.sorted.bed \
 -O z -o /scratch/rjp5nc/UK2022_2024/buscoanalysis/filteredvcf/${vcf}_busco.vcf.gz \
 /scratch/rjp5nc/UK2022_2024/daphnia_phylo/trimmed10bp_vcf/$vcf