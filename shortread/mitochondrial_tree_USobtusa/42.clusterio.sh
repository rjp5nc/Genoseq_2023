#!/usr/bin/env bash

#SBATCH -J clusterio # A single job name for the array
#SBATCH --ntasks-per-node=10 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 72:00:00 ### 15 seconds
#SBATCH --mem 50G
#SBATCH -o /scratch/rjp5nc/erroroutputs/cluster_1.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/erroroutputs/cluster_1.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab
#SBATCH --mail-type=END               # Send email at job completion
#SBATCH --mail-user=rjp5nc@virginia.edu    # Email address for notifications

#conda create -n clustalo-env -c bioconda -c conda-forge clustalo
conda activate clustalo-env

#mkdir -p /scratch/rjp5nc/UK2022_2024/consensusmitoaligned
#cat /scratch/rjp5nc/UK2022_2024/consensusmito/*.fa > /scratch/rjp5nc/UK2022_2024/consensusmitoaligned/all_mitosequences.fasta


clustalo -i /scratch/rjp5nc/UK2022_2024/consensusmitoaligned/all_mitosequences.fasta \
-o /scratch/rjp5nc/UK2022_2024/consensusmitoaligned/all_aligned.fasta --threads 10 --force
