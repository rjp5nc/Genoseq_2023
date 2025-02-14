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

### sbatch ~/Genoseq_2023/ElvisLongRead/Canu.sh
# Working directory
wd="/scratch/rjp5nc/HMW/HMWDNAElvis3"
cd ${wd}

# Modules
module load bioconda/py3.10
conda init
conda activate hicanu
# Install canu once
# conda create -n hicanu
# module load mamba
# mamba install canu



# Modules
#module load bioconda/py3.10 
#bash Miniconda3-latest-Linux-x86_64.sh
#source ~/.bashrc
#conda create --prefix ~/my_conda_env -c bioconda canu
#conda activate ~/my_conda_env

# Start
echo "Starting HiCanu"
date

# Run HiCanu
canu \
-assemble \
-p dap \
 -d dap_hifi_trim5 \
 maxThreads=40 \
 maxMemory=200g \
 useGrid=false \
 genomeSize=200m \
 gnuplot=/home/rjp5nc/miniconda3/bin/gnuplot \
 -pacbio-raw m84128_250121_222443_s2.hifi_reads.bc2104.fq.gz

### gridOptions="--partition=standard --account=berglandlab --time=2-00:00"\
### gridEngineArrayMaxJobs=5 \

# Finish
echo "Finish HiCanu"
date


