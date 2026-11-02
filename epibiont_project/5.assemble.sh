#!/usr/bin/env bash
#
#SBATCH -J megahit_full
#SBATCH --ntasks-per-node=40
#SBATCH -N 1
#SBATCH -t 72:00:00
#SBATCH --mem 200G
#SBATCH -o /scratch/rjp5nc/err/megahit_full.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/err/megahit_full.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab_standard
#SBATCH --mail-type=END               # Send email at job completion
#SBATCH --mail-user=rjp5nc@virginia.edu    # Email address for notifications

set -euo pipefail
module load miniforge
export PATH=$HOME/MEGAHIT-1.2.9-Linux-x86_64-static/bin:$PATH

cd /scratch/rjp5nc/vorticella

megahit -1 nonDaphnia_R1.fq.gz -2 nonDaphnia_R2.fq.gz \
        -o megahit_out \
        -t 40 \
        -m 0.9 \
        --min-contig-len 500