#!/usr/bin/env bash
#
#SBATCH -J basicStats # A single job name for the array
#SBATCH --ntasks-per-node=10 # ten core
#SBATCH -N 1 # on one node
#SBATCH -t 0-5:00  ### 48 hours
#SBATCH --mem 60G
#SBATCH -o /scratch/rjp5nc/basicstats.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/basicstats.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# sbatch ~/Genoseq_2023/basicStats.sh
### sacct -j 45333345


module load gcc/7.1.0  openmpi/3.1.4 R/4.1.1

Rscript --vanilla /home/rjp5nc/Genoseq_2023/basicStats.R
