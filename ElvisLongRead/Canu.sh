#!/usr/bin/env bash
#
#SBATCH -J run_hicanu # A single job name for the array
#SBATCH --ntasks-per-node=20 # multi core
#SBATCH -N 1 # on one node
#SBATCH -t 3-00:00 # 3 days
#SBATCH --mem 100G
#SBATCH -o /scratch/rjp5nc/run_hiCanu1.out # Standard output
#SBATCH -e /scratch/rjp5nc/run_hiCanu1.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# Working directory
wd="/scratch/rjp5nc/HMW/HMWDNAElvis3"
cd ${wd}

# Modules
module load bioconda/py3.10 
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc
conda create --prefix ~/my_conda_env -c bioconda canu
conda activate ~/my_conda_env

# Install canu once
# conda create -n hicanu
# module load mamba/2.0.5
# mamba install canu

# Start
echo "Starting HiCanu"
date

# Run HiCanu
canu \
-assemble \
-p dap \
 -d dap_hifi_trim \
 maxThreads=20 \
 maxMemory=100g \
 useGrid=false \
 genomeSize=150m \
 -pacbio-raw m84128_250121_222443_s2.hifi_reads.bc2104.fq.gz

# Finish
echo "Finish HiCanu"
date


