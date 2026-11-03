#!/usr/bin/env bash
#
#SBATCH -J megahit_v2
#SBATCH --ntasks-per-node=40
#SBATCH -N 1
#SBATCH -t 72:00:00
#SBATCH --mem 300G
#SBATCH -o /scratch/rjp5nc/err/megahit_v3.%A_%a.out
#SBATCH -e /scratch/rjp5nc/err/megahit_v3.%A_%a.err
#SBATCH -p standard
#SBATCH --account berglandlab_standard
#SBATCH --mail-type=END
#SBATCH --mail-user=rjp5nc@virginia.edu


set -euo pipefail
module load miniforge
export PATH=$HOME/MEGAHIT-1.2.9-Linux-x86_64-static/bin:$PATH

cd /scratch/rjp5nc/vorticella

for f in trim2_R1.fq.gz trim2_R2.fq.gz; do
  [ -s "$f" ] || { echo "Missing $f"; exit 1; }
done

megahit -1 trim2_R1.fq.gz -2 trim2_R2.fq.gz \
        -o megahit_out3 \
        -t 40 \
        -m 0.9 \
        --min-contig-len 500