#!/usr/bin/env bash
#
#SBATCH -J RemoveRepeats # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-2:00  ### 1 hours
#SBATCH --mem 24G
#SBATCH -o /scratch/rjp5nc/UK2022_2024/allshortreads/chr/2022seq.concat.Removereps.renamed.vcf.gz
#SBATCH -e /scratch/rjp5nc/outputerrors/download_%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

### sbatch ~/Genoseq_2023/RemoveRepetitiveRegions.sh
### sacct -j 22867938
### scancel
### cat

module load gcc/12.4.0 bedtools/2.30.0 htslib

bedtools subtract -header \
-a /scratch/rjp5nc/UK2022_2024/allshortreads/chr/2022seq.concat.renamed.vcf.gz \
-b /project/berglandlab/daphnia_ref/RMoutHiCGMgoodscaff.keep.bed