#!/usr/bin/env bash
#
#SBATCH -J touch # A single job name for the array
#SBATCH --ntasks-per-node=2 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-48:00  ### 48 hours
#SBATCH --mem 10G
#SBATCH -o /scratch/rjp5nc/err/move.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/err/move.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

cd /scratch/rjp5nc
find . -type f -exec touch {} +
find . -type d -exec touch {} +
