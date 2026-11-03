#!/usr/bin/env bash
#
#SBATCH -J vort_persample
#SBATCH --ntasks-per-node=16
#SBATCH -N 1
#SBATCH -t 12:00:00
#SBATCH --mem 32G
#SBATCH --array=1-36
#SBATCH -o /scratch/rjp5nc/err/vort_persample.%A_%a.out
#SBATCH -e /scratch/rjp5nc/err/vort_persample.%A_%a.err
#SBATCH -p standard
#SBATCH --account berglandlab_standard

set -euo pipefail
module purge
module load bowtie2 samtools

raw=/scratch/rjp5nc/rawdata/mysamples
wd=/scratch/rjp5nc/vorticella
idx=$wd/contigs_v2_idx/contigs_v2
out=$wd/persample
mkdir -p "$out"

# Sample for this array task: Rockpool4 D1..F12
rows=(D E F)
i=$((SLURM_ARRAY_TASK_ID - 1))
s=Rockpool4_${rows[$((i / 12))]}$(( i % 12 + 1 ))

r1=$(ls "$raw"/"$s"/*_1.fq.gz | paste -sd,)
r2=$(ls "$raw"/"$s"/*_2.fq.gz | paste -sd,)

bowtie2 -x "$idx" -1 "$r1" -2 "$r2" -p 16 --local --no-unal 2> "$out"/"$s".log \
  | samtools sort -@ 4 -o "$out"/"$s".bam -
samtools index "$out"/"$s".bam
samtools idxstats "$out"/"$s".bam > "$out"/"$s".idxstats

# Total input reads (from bowtie2 log) and reads on Vorticella bin contigs
total=$(head -1 "$out"/"$s".log | awk '{print $1*2}')
vort=$(awk 'NR==FNR{k[$1]=1; next} ($1 in k){s+=$3} END{print s+0}' \
  "$wd"/vort_bin/vort_bin2_ids.txt "$out"/"$s".idxstats)
echo -e "$s\t$total\t$vort" > "$out"/"$s".vortcount