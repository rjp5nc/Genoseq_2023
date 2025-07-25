#!/usr/bin/env bash

#SBATCH -J clustalo # A single job name for the array
#SBATCH --ntasks-per-node=10 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 72:00:00 ### 15 seconds
#SBATCH --mem 20G
#SBATCH -o /scratch/rjp5nc/erroroutputs/cluster_1.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/erroroutputs/cluster_1.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab
#SBATCH --mail-type=END               # Send email at job completion
#SBATCH --mail-user=rjp5nc@virginia.edu    # Email address for notifications

#conda create -n clustalo-env -c bioconda -c conda-forge clustalo
source ~/miniconda3/etc/profile.d/conda.sh
conda activate clustalo-env
#conda install -c bioconda fasttree

#mkdir -p /scratch/rjp5nc/UK2022_2024/consensusmitoaligned
cat /scratch/rjp5nc/UK2022_2024/consensusmito/*.fa > /scratch/rjp5nc/UK2022_2024/consensusmitoaligned/all_mitosequences.fasta

clustalo -i /scratch/rjp5nc/UK2022_2024/consensusmitoaligned/all_mitosequences.fasta \
-o /scratch/rjp5nc/UK2022_2024/consensusmitoaligned/all_aligned.fasta --threads 10 --force

awk '/^>/{count[$0]++; if(count[$0]>1) $0=$0"_"count[$0]; print; next}1' \
/scratch/rjp5nc/UK2022_2024/consensusmitoaligned/all_aligned.fasta \
> /scratch/rjp5nc/UK2022_2024/consensusmitoaligned/all_aligned_unique.fasta

FastTree -nt -gtr /scratch/rjp5nc/UK2022_2024/consensusmitoaligned/all_aligned_unique.fasta > /scratch/rjp5nc/UK2022_2024/consensusmitoaligned/mito_tree.nwk

#awk 'BEGIN{n=0} /^>/ {n++} n<=100 {print}' /scratch/rjp5nc/UK2022_2024/consensusmitoaligned/all_mitosequences.fasta > mito_subset.fasta

# Run locally with 4 threads (or adjust to 10)
#clustalo -i mito_subset.fasta -o mito_subset.aln --threads=4