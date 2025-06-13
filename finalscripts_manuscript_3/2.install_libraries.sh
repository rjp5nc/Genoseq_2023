#!/usr/bin/env bash

#SBATCH -J runFASTQC # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 48:00:00 ### 15 seconds
#SBATCH --mem 6G
#SBATCH -o /scratch/rjp5nc/Canu_error/nFlo_1.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/Canu_error/nFlo_1.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# sbatch /home/rjp5nc/Genoseq_2023/ElvisLongRead/install_libraries.sh
# sacct -j 2846857
cd /scratch/rjp5nc/krakenDB
wget https://genome-idx.s3.amazonaws.com/kraken/k2_core_nt_20241228.tar.gz
wget https://genome-idx.s3.amazonaws.com/kraken/core_nt_20241228/inspect.txt
wget https://genome-idx.s3.amazonaws.com/kraken/core_nt_20241228/library_report.tsv