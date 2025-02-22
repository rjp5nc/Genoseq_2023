#!/usr/bin/env bash
#
#SBATCH -J run_hicanu2 # A single job name for the array
#SBATCH --ntasks-per-node=40 # multi core
#SBATCH -N 1 # on one node
#SBATCH -t 6-00:00 # 5 days
#SBATCH --mem 200G
#SBATCH -o /scratch/rjp5nc/Canu_error/run_hiCanu1%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/Canu_error/run_hiCanu1%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

### sbatch ~/Genoseq_2023/ElvisLongRead/Canu2.sh
# Working directory
wd="/scratch/rjp5nc/assemblecontigs"
cd ${wd}

# Modules
module load bioconda/py3.10 samtools
source ~/.bashrc
conda activate hicanu

# Install canu once
# conda create -n hicanu
# module load mamba
# mamba install canu


# Start
echo "Starting HiCanu"
date

# Run HiCanu
canu \
-p asm \
 -d asm_output \
 maxThreads=40 \
 maxMemory=200g \
 useGrid=false \
 genomeSize=200m \
 gnuplot=/home/rjp5nc/miniconda3/bin/gnuplot \
 -pacbio-hifi dap.contigs.fasta

# Finish
echo "Finish HiCanu"
date


