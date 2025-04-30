#!/usr/bin/env bash

#SBATCH -J runFASTQC # A single job name for the array
#SBATCH --ntasks-per-node=20 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 4:00:00 ### 15 seconds
#SBATCH --mem 120G
#SBATCH -o /scratch/rjp5nc/erroroutputs/dict.out # Standard output
#SBATCH -e /scratch/rjp5nc/erroroutputs/dict.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

module load R/4.4.2;R

Rscript --vanilla /home/rjp5nc/Genoseq_2023/shortread/liftover/masterdictionary.R
