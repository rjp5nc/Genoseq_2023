#!/usr/bin/env bash
#
#SBATCH -J vcf2gds # A single job name for the array
#SBATCH --ntasks-per-node=10 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-48:00  ### 48 hours
#SBATCH --mem 24G
#SBATCH -o /scratch/rjp5nc/vcf2gds.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/vcf2gds.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

### sbatch ~/Dappupool20182019/variantCalling/vcf2gds.sh
### sacct -j 22867938
### cat /scratch/aob2x/dest/slurmOutput/vcf2gds.22867938

module load gcc/7.1.0 openmpi/3.1.4 R/4.1.1

Rscript --vanilla /home/rjp5nc/Dappupool20182019/variantCalling/vcf2gds.R \
/project/berglandlab/Robert/DorsetPooledSequencing2018_2019/AllPooledVCFsIncludingDcatDbunk/Allpools.concat.Removereps.vcf.gz
