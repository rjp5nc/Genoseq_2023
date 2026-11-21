#!/usr/bin/env bash
#
#SBATCH -J phyloflash_A1
#SBATCH --ntasks-per-node=16
#SBATCH -N 1
#SBATCH -t 12:00:00
#SBATCH --mem 48G
#SBATCH -o /scratch/rjp5nc/err/phyloflash_A1.%A_%a.out
#SBATCH -e /scratch/rjp5nc/err/phyloflash_A1.%A_%a.err
#SBATCH -p standard
#SBATCH --account berglandlab_standard

set -euo pipefail
export PATH=/home/rjp5nc/miniconda3/envs/phyloflash/bin:$PATH

sample=Gilmer5_A1
raw=/scratch/rjp5nc/rawdata/mysamples/$sample
out=/scratch/rjp5nc/ssu_analysis/phyloflash/$sample
db=/scratch/rjp5nc/phyloflash_db/138.1.eukmod
mkdir -p "$out"
cd "$out"

# Combine all lanes (same order in R1 and R2 keeps pairs matched)
cat "$raw"/*_1.fq.gz > ${sample}_R1.fq.gz
cat "$raw"/*_2.fq.gz > ${sample}_R2.fq.gz

phyloFlash.pl -lib "$sample" \
  -read1 ${sample}_R1.fq.gz -read2 ${sample}_R2.fq.gz \
  -readlength 150 -CPUs 16 \
  -dbhome "$db" \
  -taxlevel 6 \
  -almosteverything

rm -f ${sample}_R1.fq.gz ${sample}_R2.fq.gz
echo "Done: $(date)"