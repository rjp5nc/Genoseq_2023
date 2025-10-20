#!/usr/bin/env bash

#SBATCH -J BEAST # A single job name for the array
#SBATCH --ntasks-per-node=30 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 6-0:00:00 ### 15 seconds
#SBATCH --mem 200G
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

cd /scratch/rjp5nc/snapp5/

#/scratch/rjp5nc/beast/beast/bin/packagemanager -add SNAPP

/scratch/rjp5nc/beast/beast/bin/beast -java -seed 12354 -threads 25 /scratch/rjp5nc/snapp5/snapp.mito_only1.xml

#Ran treeannotator locally with increased RAM.


