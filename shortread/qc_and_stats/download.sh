#!/usr/bin/env bash
#
#SBATCH -J download2l # A single job name for the array
#SBATCH --ntasks-per-node=10 # ten core
#SBATCH -N 1 # on one node
#SBATCH -t 0-36:00  ### 24 hours
#SBATCH --mem 10G
#SBATCH -o /scratch/rjp5nc/basicstats.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/basicstats.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# sbatch ~/Genoseq_2023/download.sh
### sacct -j 45333345 



module load lftp

cd /scratch/rjp5nc/2ldup

lftp -e "set sftp:auto-confirm yes; open sftp://X202SC25016886-Z01-F001:dxk06fch@usftp23.novogene.com:3022; mirror --verbose --use-pget-n=8 -c"