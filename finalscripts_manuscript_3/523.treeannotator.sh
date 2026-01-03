#!/usr/bin/env bash

#SBATCH -J treeannotator # A single job name for the array
#SBATCH --ntasks-per-node=20 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 1:00:00 ### 15 seconds
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

cd /scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/gatk_gvcf_top3_plusA/

ls /scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/gatk_gvcf_top3_plusA/

# Open resutls with Tracer
#/home/rjp5nc/tracer/bin/tracer

# Open density tree with beast2
#/home/rjp5nc/beast/bin/densitree

/scratch/rjp5nc/beast/beast/bin/treeannotator \
-burnin 10 \
-height mean \
-topology MCC \
/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/gatk_gvcf_top3/snapp.mito.top3_plusA.trees \
/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/gatk_gvcf_top3_plusA/snapp.mito.highest1.tree

