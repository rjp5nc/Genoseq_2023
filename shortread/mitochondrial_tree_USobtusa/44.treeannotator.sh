#!/usr/bin/env bash

#SBATCH -J treeannotator # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 6-0:00:00 ### 15 seconds
#SBATCH --mem 85G
#SBATCH -o /scratch/rjp5nc/erroroutputs/beast.%A.out # Standard output
#SBATCH -e /scratch/rjp5nc/erroroutputs/beast.%A.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab
#SBATCH --mail-type=END               # Send email at job completion
#SBATCH --mail-user=rjp5nc@virginia.edu    # Email address for notifications

#conda create -n beast2-277 openjdk=17 -c conda-forge -y
source ~/miniconda3/etc/profile.d/conda.sh
conda activate beast2-277

cd /scratch/rjp5nc/UK2022_2024/consensusmitoaligned2/

java -Xmx70G -cp "/scratch/rjp5nc/beast/beast/lib/*" beast.app.tools.TreeAnnotator \
  -burnin 10 \
  -heights mean \
  all_aligned-all_aligned_unique.trees \
  annotated_MCC.tree

#beastfx.app.treeannotator.TreeAnnotator

