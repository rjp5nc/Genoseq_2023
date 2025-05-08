#!/usr/bin/env bash
#
#SBATCH -J vcf2gds
#SBATCH --ntasks-per-node=10 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 1-00:00 # hours
#SBATCH --mem 60G
#SBATCH -o /scratch/rjp5nc/err/gatk.chrom.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/err/gatk.chrom.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# Modules to load
module load goolf/11.4.0_4.1.4 R/4.5.0; module load gdal geos proj

# Start
echo date

# Run script
Rscript /home/rjp5nc/Genoseq_2023/shortread/whole_dataset/vcf2gds.R

# Finish
echo "Finish"
echo date