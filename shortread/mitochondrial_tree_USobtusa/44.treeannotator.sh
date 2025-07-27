#!/usr/bin/env bash

#SBATCH -J treeannotator # A single job name for the array
#SBATCH --ntasks-per-node=30 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 6-0:00:00 ### 15 seconds
#SBATCH --mem 120G
#SBATCH -o /scratch/rjp5nc/erroroutputs/beast.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/erroroutputs/beast.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab
#SBATCH --mail-type=END               # Send email at job completion
#SBATCH --mail-user=rjp5nc@virginia.edu    # Email address for notifications

#conda create -n beast2-277 openjdk=17 -c conda-forge -y
source ~/miniconda3/etc/profile.d/conda.sh
conda activate beast2-277

cd /scratch/rjp5nc/UK2022_2024/consensusmitoaligned2/

/scratch/rjp5nc/beast/beast/bin/treeannotator \
  -burnin 10 \
  -heights mean \
  all_aligned-all_aligned_unique.trees \
  annotated_MCC.tree
