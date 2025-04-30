#!/usr/bin/env bash

#SBATCH -J crossmap # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 12:00:00 
#SBATCH --mem 200G
#SBATCH -o /scratch/rjp5nc/erroroutputs/crossmap.out # Standard output
#SBATCH -e /scratch/rjp5nc/erroroutputs/crossmap.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab
#SBATCH --mail-type=END               # Send email at job completion
#SBATCH --mail-user=rjp5nc@virginia.edu    # Email address for notifications

#Load necessary modules (if needed)
module load gcc/11.4.0
module load openmpi/4.1.4
module load python/3.11.4

#pip install --user .

#species=us_obtusa
#species=us_ambigua
#species=eu_obtusa
species=us_pulex

bed=us_pulex_ref_kap4.allbases.cleaned.slimmed.bed

/home/rjp5nc/.local/bin/CrossMap bed --chromid a \
/scratch/rjp5nc/lastz/$species/chainnet/liftover.chain \
/scratch/rjp5nc/Reference_genomes/$bed \
/scratch/rjp5nc/liftover/uspulex_to_eupulex.bed

#cd /home/rjp5nc/Genoseq_2023/shortread/liftover