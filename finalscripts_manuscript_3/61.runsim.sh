#!/usr/bin/env bash
#
#SBATCH -J stats # A single job name for the array
#SBATCH --ntasks-per-node=10 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-72:00  ### 48 hours
#SBATCH --mem 60G
#SBATCH -o /scratch/rjp5nc/err/stats.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/err/stats.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab


module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1;R

Rscript --vanilla /home/rjp5nc/Genoseq_2023/finalscripts_manuscript_3/61.genomic.similarity_fst.R