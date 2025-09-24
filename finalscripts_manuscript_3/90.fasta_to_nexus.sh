#!/usr/bin/env bash
#
#SBATCH -J trim # A single job name for the array
#SBATCH --ntasks-per-node=10 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-10:00  ### 48 hours
#SBATCH --mem 50G
#SBATCH -o /scratch/rjp5nc/err/trim.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/err/trim.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab


cd /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/

cat consensus_fastas/*.fasta > all_samples.fasta
python /scratch/rjp5nc/AMAS/amas/AMAS.py convert -i all_samples.fasta -f fasta -d dna -o all_samples.nex
