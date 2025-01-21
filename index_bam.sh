#!/usr/bin/env bash
#
#SBATCH -J download # A single job name for the array
#SBATCH --ntasks-per-node=10 # ten cores
#SBATCH -N 1 # on one node
#SBATCH -t 1:00:00 ### 8 hours
#SBATCH --mem 60G
#SBATCH -o /scratch/rjp5nc/download_%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/download_%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

module load samtools/1.17

# sbatch --array=1-23 ~/Dappupool20182019/Mapping/index_bam.sh
### sacct -j 40641531
###  cat /scratch/rjp5nc/download_40446865_1.err
# sacct -j

### example:
# SLURM_ARRAY_TASK_ID=1

dir=$1
echo ${dir}

samtools index \
-@ 10 \
${dir} \
${dir}.bai
