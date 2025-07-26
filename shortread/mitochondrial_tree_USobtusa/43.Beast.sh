#!/usr/bin/env bash

#SBATCH -J BEAST # A single job name for the array
#SBATCH --ntasks-per-node=30 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 6-0:00:00 ### 15 seconds
#SBATCH --mem 120G
#SBATCH -o /scratch/rjp5nc/erroroutputs/beast.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/erroroutputs/beast.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab
#SBATCH --mail-type=END               # Send email at job completion
#SBATCH --mail-user=rjp5nc@virginia.edu    # Email address for notifications

#conda create -n beast2-277 openjdk=17 -c conda-forge -y
source ~/miniconda3/etc/profile.d/conda.sh
conda activate beast2-277

#wget https://github.com/CompEvol/beast2/releases/download/v2.7.7/BEAST.v2.7.7.Linux.tgz
#tar -xzf BEAST.v2.7.7.Linux.tgz

#conda install -c bioconda libbeagle

cd /scratch/rjp5nc/UK2022_2024/consensusmitoaligned2/

#I had to increase the ram usage in the .bat file of Beauti2 to turn the alignment into a .xml file. It wasnt saving otherwise

#beast -beagle_SSE -threads 30 /scratch/rjp5nc/UK2022_2024/consensusmitoaligned2/all_aligned.xml

beast -beagle_SSE -beagle_instances 30 -threads 30 /scratch/rjp5nc/UK2022_2024/consensusmitoaligned2/all_aligned.xml
