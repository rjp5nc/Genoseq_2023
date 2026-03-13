#!/usr/bin/env bash
#
#SBATCH -J stats # A single job name for the array
#SBATCH --ntasks-per-node=10 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 4-00:00  ### 48 hours
#SBATCH --mem 150G
#SBATCH -o /scratch/rjp5nc/err/stats.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/err/stats.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

### sbatch ~/Genoseq_2023/shortread/liftover/depth_GDS.sh
### sacct -j 22867938

module load gcc openmpi R/4.3.1

Rscript --vanilla /home/rjp5nc/Genoseq_2023/finalscripts_manuscript_3/45.stats_busco+1.R