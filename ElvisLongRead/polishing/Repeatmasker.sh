#!/usr/bin/env bash
#
#SBATCH -J run_masker # A single job name for the array
#SBATCH --ntasks-per-node=10 # multi core
#SBATCH -N 1 # on one node
#SBATCH -t 1-00:00 # 1 days
#SBATCH --mem 150G
#SBATCH -o /scratch/rjp5nc/err/m%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/err/m%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab


conda activate repeatmodeler_new

path=/scratch/rjp5nc/removedups/finalfiles/usobtusa

mkdir -p $path

cd $path 

ref=/scratch/rjp5nc/Reference_genomes/post_kraken/US_obtusa_onlydaps.fa
classified=/scratch/rjp5nc/removedups/us_dobtusa/RM_653473.SunMar301757312025/consensi.fa.classified

export REPEATMASKER_NO_GENLIB=1
export REPEATMASKER_NO_TAXONOMY=1

RepeatMasker -lib $classified -dir . $ref