#!/usr/bin/env bash
#
#SBATCH -J move # A single job name for the array
#SBATCH --ntasks-per-node=10 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-48:00  ### 48 hours
#SBATCH --mem 60G
#SBATCH -o /scratch/rjp5nc/err/move.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/err/move.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab


mkdir -p /scratch/rjp5nc/rawdata
cp -r /project/berglandlab/Robert/shortread_data/data/01.RawData/* /scratch/rjp5nc/rawdata/
