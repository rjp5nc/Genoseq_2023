#!/usr/bin/env bash
#
#SBATCH -J move # A single job name for the array
#SBATCH --ntasks-per-node=40 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 72:00:00
#SBATCH --mem 200G
#SBATCH -o /scratch/rjp5nc/err/move.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/err/move.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab_standard
#SBATCH --mail-type=END               # Send email at job completion
#SBATCH --mail-user=rjp5nc@virginia.edu    # Email address for notifications

set -euo pipefail
export PATH=$HOME/MEGAHIT-1.2.9-Linux-x86_64-static/bin:$PATH

module load seqkit

cd /scratch/rjp5nc/vorticella

# Random 10% of pairs; same seed on both files keeps pairs matched
seqkit sample -p 0.1 -s 42 -j 8 nonDaphnia_R1.fq.gz -o sub10_R1.fq.gz
seqkit sample -p 0.1 -s 42 -j 8 nonDaphnia_R2.fq.gz -o sub10_R2.fq.gz
seqkit stats sub10_R1.fq.gz sub10_R2.fq.gz

megahit -1 sub10_R1.fq.gz -2 sub10_R2.fq.gz \
        -o megahit_sub10 \
        -t 40 -m 0.9 \
        --min-contig-len 500



# cd /scratch/rjp5nc/vorticella

# megahit -1 nonDaphnia_R1.fq.gz -2 nonDaphnia_R2.fq.gz \
#         -o megahit_out \
#         -t 40 \
#         -m 0.9 \
#         --min-contig-len 500