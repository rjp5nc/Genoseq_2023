#!/usr/bin/env bash
#
#SBATCH -J run_earl_grey # A single job name for the array
#SBATCH --ntasks-per-node=10 # multi core
#SBATCH -N 1 # on one node
#SBATCH -t 1-00:00 # 1 days
#SBATCH --mem 40G
#SBATCH -o /scratch/rjp5nc/err/m%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/err/m%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

source ~/miniconda3/etc/profile.d/conda.sh

#conda create -n earlgrey_env python=3.9
conda activate earlgrey_env

#conda install -c bioconda -c conda-forge repeatmasker=4.1.8
#CONDA_VERBOSE=3 conda install -c bioconda -c conda-forge earlgrey

cd /scratch/rjp5nc/UK2022_2024/daphnia_phylo/earlgrey

ref=/scratch/rjp5nc/Reference_genomes/post_kraken/assembly.hap2_onlydaps.fasta
classified=/scratch/rjp5nc/removedups/eu_dobtusa/RM_935821.MonMar240837102025/consensi.fa.classified

/home/rjp5nc/miniconda3/envs/earlgrey_env/share/earlgrey-4.2.4-0/earlGrey -g $ref -s Daphniaobtusa -o earlgrey_output -l $classified

#consensi.fa.classified from repeatmodeller
#genome.fasta is the assembled fasta

