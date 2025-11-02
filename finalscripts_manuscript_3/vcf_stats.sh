#!/usr/bin/env bash
#
#SBATCH -J vcfstats
#SBATCH --ntasks-per-node=4 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-72:00 # hours
#SBATCH --mem 40G
#SBATCH -o /scratch/rjp5nc/err/FilterVCFs.out # Standard output
#SBATCH -e /scratch/rjp5nc/err/FilterVCFs.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

module load bcftools/1.17

VCF="lifted_12major.vcf.gz"
out="lifted_12major_stats.vchk"

cd /scratch/rjp5nc/UK2022_2024/daphnia_phylo/trimmed_10bp_repeatmasked_vcf/lifted_vcf

bcftools stats $VCF > $out

mkdir -p plots

#conda create -n vcfplots python=3.9 matplotlib -y
conda activate vcfplots
#cd plots/
#python plot.py

plot-vcfstats -p plots/ $out
