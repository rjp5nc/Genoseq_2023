#!/usr/bin/env bash
#
#SBATCH -J move # A single job name for the array
#SBATCH --ntasks-per-node=10 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-48:00  ### 48 hours
#SBATCH --mem 60G
#SBATCH -o /scratch/rjp5nc/move.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/move.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab


##find /project/berglandlab/connor -type f -name "genomedb"
##/project/berglandlab/connor/backup_project/BACKUP_scratch/DBI/

##rsync -a --prune-empty-dirs --include "*/" --include "*_1.fq.gz" --include "*_2.fq.gz" --exclude "*" /scratch/rjp5nc/UK2022_2024/allshortreads/01.RawData /project/berglandlab/Robert/shortread_data/data


cp /project/berglandlab/connor/backup_project/BACKUP_scratch/DBI /scratch/rjp5nc