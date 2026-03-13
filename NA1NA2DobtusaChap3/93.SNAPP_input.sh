#!/usr/bin/env bash
#
#SBATCH -J Snapp_input # A single job name for the array
#SBATCH --ntasks-per-node=15 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 6-0:00  ### 48 hours
#SBATCH --mem 200G
#SBATCH -o /scratch/rjp5nc/err/trim.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/err/trim.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

#!/bin/bash
# Script used to prepare VCF and run snapp for phylogenetic splits

# Load modules
module load bcftools
module load ruby

# Working directory
wd="/scratch/rjp5nc/snapp5"

cd $wd

# VCF
vcf="/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/trimmed10bp_filtered_two_of_each_fixed2.vcf.gz"

# Beast2 directory
beast2="/scratch/rjp5nc/beast/beast/bin/beast"


IN=/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/two_of_each_genomic_type.csv
OUTTMP=/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/two_of_each_clone_genomic.tmp.txt
OUT=/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/two_of_each_clone_genomic.txt
OUTOO=/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/two_of_each_clone_genomic_unique_OO.txt

awk -F, 'NR==1 {print "CloneA\tGenomic_type"; next} {gsub(/"/,""); print $4"\t"$2}' $IN > $OUTTMP
tail -n +2 $OUTTMP > $OUT

awk '{
    if($2=="OO"){ 
        count++; 
        $2="OO_"count 
    } 
    print $0
}' $OUT > $OUTOO


awk '{print $2, $1}' /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/two_of_each_clone_genomic_unique_OO.txt > /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/two_of_each_clone_genomic_unique_OO_swapped.txt

# 5 indviduals per species - highest mean coverage
samps=/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/final_vcf_filter_two_of_each.txt

# Starting tree: base_tree.nwk
#(((((Daphnia.pulex.NorthAmerica,Daphnia.pulexcaria.NorthAmerica,Daphnia.pulicaria.NorthAmerica),Daphnia.pulicaria.Europe),Daphnia.pulex.Europe),Daphnia.obtusa.NorthAmerica),Daphnia.obtusa.Europe)

# Only biallelic SNPs 
#bcftools view \
#--threads 15 \
#-e 'AC==0 || AC==AN' \
#-m2 \
#-M2 \
#-Ov \
#-o ${wd}/daphnia.genome.2inds.biallelic.vcf \
#${vcf}

#bcftools view \
#  -e 'ALT="*"' \
#  -m2 -M2 \
#  -Ov \
#  -o ${wd}/daphnia.genome.2inds.biallelic.clean.vcf \
#  ${vcf}

#cat /project/berglandlab/connor/chapter2_TSP/divergence_snps/snapp5/constraints.withobtusa.txt

# Run ruby input script
# From: https://github.com/mmatschiner/snapp_prep
ruby ${wd}/snapp_prep.rb \
-v ${wd}/daphnia.genome.2inds.biallelic.clean.vcf \
-t /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/two_of_each_clone_genomic_unique_OO_swapped.txt \
-c /scratch/rjp5nc/snapp5/constraints.obtusa_monomorphic.txt \
-x ${wd}/snapp.mono2.xml \
-o ${wd}/snapp.mono2 \
-m 1000 \
-l 100000

# Open resutls with Tracer
#/home/csm6hg/tracer/bin/tracer

# Open density tree with beast2
#/home/csm6hg/beast/bin/densitree