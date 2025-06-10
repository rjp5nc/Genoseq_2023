#!/usr/bin/env bash
#
#SBATCH -J run_masker # A single job name for the array
#SBATCH --ntasks-per-node=10 # multi core
#SBATCH -N 1 # on one node
#SBATCH -t 1-00:00 # 1 days
#SBATCH --mem 150G
#SBATCH -o /scratch/rjp5nc/err/m%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/err/m%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

source ~/miniconda3/etc/profile.d/conda.sh
conda activate repeatmodeler_new

#us pulex, us ambigua, eu obtusa, eu pulex

path=/scratch/rjp5nc/removedups/finalfiles/usambigua

mkdir -p $path

cd $path 

##US_obtusa_onlydaps.fa
##us_pulex_ref_kap4.fa
#assembly.hap2_onlydaps.fasta
##Daphnia_ambigua_Q001_genome.fa
#totalHiCwithallbestgapclosed.fa

ref=/scratch/rjp5nc/Reference_genomes/post_kraken/Daphnia_ambigua_Q001_genome.fa
classified=/scratch/rjp5nc/removedups/us_dambigua/RM_229818.SunMar301804372025/consensi.fa.classified

RepeatMasker -lib $classified -dir . $ref

#cd /home/rjp5nc/Genoseq_2023/ElvisLongRead/polishing