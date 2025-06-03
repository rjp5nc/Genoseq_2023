#!/usr/bin/env bash
#
#SBATCH -J depth # A single job name for the array
#SBATCH --ntasks-per-node=10 # ten core
#SBATCH -N 1 # on one node
#SBATCH -t 0-10:00  ### 48 hours
#SBATCH --mem 100G
#SBATCH -o /scratch/rjp5nc/basicstats.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/basicstats.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# sbatch ~/Genoseq_2023/tree_building.sh
### sacct -j 45333345


module load gcc/11.4.0  openmpi/4.1.4 icu R/4.3.1

Rscript --vanilla /home/rjp5nc/Genoseq_2023/ElvisLongRead/polishing/depth.r
