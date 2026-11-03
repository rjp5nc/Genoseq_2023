#!/usr/bin/env bash
#
#SBATCH -J mitofinder_megahit
#SBATCH --ntasks-per-node=40
#SBATCH -N 1
#SBATCH -t 48:00:00
#SBATCH --mem 100G
#SBATCH -o /scratch/rjp5nc/err/mitofinder_megahit.%A_%a.out
#SBATCH -e /scratch/rjp5nc/err/mitofinder_megahit.%A_%a.err
#SBATCH -p standard
#SBATCH --account berglandlab_standard
#SBATCH --mail-type=END
#SBATCH --mail-user=rjp5nc@virginia.edu

set -euo pipefail
cd /scratch/rjp5nc/vorticella/mitofinder


export PATH=/home/rjp5nc/miniconda3/envs/mitofinder/bin:$PATH
mitofinder --version

mitofinder -j vort_mito_megahit \
    -a megahit_all.fa \
    -r ciliate_mito_refs_clean.gb \
    -o 4 \
    -p 40 -m 100 \
    --new-genes \
    --tRNA-annotation trnascan \
    --allow-intron

echo "Done: $(date)"