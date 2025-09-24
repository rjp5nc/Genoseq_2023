#!/usr/bin/env bash
#
#SBATCH -J trim # A single job name for the array
#SBATCH --ntasks-per-node=30 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 6-0:00  ### 48 hours
#SBATCH --mem 200G
#SBATCH -o /scratch/rjp5nc/err/trim.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/err/trim.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

module load gcc/11.4.0 openmpi/4.1.4 mafft/7.505

cd /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/

#awk '/^>/{count[$0]++; if(count[$0]>1) $0=$0"_"count[$0]; print; next}1' \
#only2each.fasta \
#> only2each_unique_notaligned.fasta

mafft --thread 30 --retree 1 --maxiterate 0 only2each_unique_notaligned.fasta > all_samples_mafft.fasta
