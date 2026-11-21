#!/usr/bin/env bash
#
#SBATCH -J phyloflash
#SBATCH --ntasks-per-node=16
#SBATCH -N 1
#SBATCH -t 6:00:00
#SBATCH --mem 48G
#SBATCH -o /scratch/rjp5nc/err/phyloflash.%A_%a.out
#SBATCH -e /scratch/rjp5nc/err/phyloflash.%A_%a.err
#SBATCH -p standard
#SBATCH --account berglandlab_standard
#SBATCH --array=1-576%40

set -euo pipefail
export PATH=/home/rjp5nc/miniconda3/envs/phyloflash/bin:$PATH

base=/scratch/rjp5nc/ssu_analysis/phyloflash
sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$base/samples.txt")
raw=/scratch/rjp5nc/rawdata/mysamples/$sample
out=$base/$sample
db=/scratch/rjp5nc/phyloflash_db/138.1.eukmod

# Skip samples that already finished
if [[ -s "$out/$sample.phyloFlash.tar.gz" ]]; then
  echo "Already done: $sample"; exit 0
fi

rm -rf "$out"; mkdir -p "$out"; cd "$out"

# Combine all lanes (identical sort order keeps R1/R2 paired)
ls "$raw"/*_1.fq.gz > /dev/null
cat $(ls "$raw"/*_1.fq.gz | sort) > ${sample}_R1.fq.gz
cat $(ls "$raw"/*_2.fq.gz | sort) > ${sample}_R2.fq.gz

phyloFlash.pl -lib "$sample" \
  -read1 ${sample}_R1.fq.gz -read2 ${sample}_R2.fq.gz \
  -readlength 150 -CPUs 16 \
  -dbhome "$db" \
  -taxlevel 6 \
  -almosteverything

rm -f ${sample}_R1.fq.gz ${sample}_R2.fq.gz
echo "Done: $sample $(date)"