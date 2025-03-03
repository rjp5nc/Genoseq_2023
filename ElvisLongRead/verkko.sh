#!/usr/bin/env bash
#
#SBATCH -J run_verkko # A single job name for the array
#SBATCH --ntasks-per-node=20 # multi core
#SBATCH -N 1 # on one node
#SBATCH -t 5-00:00 # 5 days
#SBATCH --mem 100G
#SBATCH -o /scratch/rjp5nc/Canu_error/verkko%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/Canu_error/verkko%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

#conda create -n verkko -c conda-forge -c bioconda -c defaults verkko

conda activate verkko
verkko -d /scratch/rjp5nc/HMW/HMWDNAElvis3/verkko --hifi /scratch/rjp5nc/HMW/HMWDNAElvis3/verkko/m84128_250121_222443_s2.hifi_reads.bc2104.fq.gz