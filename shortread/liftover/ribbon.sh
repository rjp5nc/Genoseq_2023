#!/usr/bin/env bash

#SBATCH -J SYNY # A single job name for the array
#SBATCH --ntasks-per-node=40 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 4:00:00 ### 15 seconds
#SBATCH --mem 120G
#SBATCH -o /scratch/rjp5nc/erroroutputs/nFlo_1.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/erroroutputs/nFlo_1.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

#conda create -n syny  

# Activating the conda environment
conda activate syny

## Installing syny within the conda environment
#conda install syny -c conda-forge -c bioconda 

cd /scratch/rjp5nc/Reference_genomes/post_kraken

#fasta_to_gbff.pl \
#  --fasta Daphnia_ambigua_Q001_genome.fa \
#  --outdir GBFF \
#  --gzip

cd GBFF

  run_syny.pl \
  -a *.gbff.gz \
  -o SYNY