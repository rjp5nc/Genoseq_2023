#!/usr/bin/env bash
#
#SBATCH -J map_contigs_v2
#SBATCH --ntasks-per-node=40
#SBATCH -N 1
#SBATCH -t 24:00:00
#SBATCH --mem 60G
#SBATCH -o /scratch/rjp5nc/err/map_contigs_v2.%A_%a.out
#SBATCH -e /scratch/rjp5nc/err/map_contigs_v2.%A_%a.err
#SBATCH -p standard
#SBATCH --account berglandlab_standard
#SBATCH --mail-type=END
#SBATCH --mail-user=rjp5nc@virginia.edu

set -euo pipefail
module load bowtie2
module load samtools

threads=40
wd=/scratch/rjp5nc/vorticella
asm=$wd/megahit_out2/final.contigs.fa
idxdir=$wd/contigs_v2_idx
idx=$idxdir/contigs_v2
out=$wd/coverage_v2

mkdir -p "$idxdir" "$out"
cd "$wd"

# 1. Build index of the assembly (skipped if it exists)
if [ ! -e "${idx}.1.bt2" ] && [ ! -e "${idx}.1.bt2l" ]; then
  bowtie2-build --threads "$threads" "$asm" "$idx"
fi

# 2. Map, keep only aligned reads, coordinate-sort
bowtie2 -x "$idx" -1 nonDaphnia2_R1.fq.gz -2 nonDaphnia2_R2.fq.gz \
        -p "$threads" --no-unal 2> "$out"/mapping.log \
  | samtools sort -@ 8 -m 2G -T "$out"/sorttmp -o "$out"/contigs_v2.sorted.bam -

# 3. Index BAM
samtools index -@ 8 "$out"/contigs_v2.sorted.bam

# 4. Summaries
samtools flagstat -@ 8 "$out"/contigs_v2.sorted.bam > "$out"/flagstat.txt
samtools coverage "$out"/contigs_v2.sorted.bam > "$out"/per_contig_coverage.tsv

cat "$out"/mapping.log


cd /scratch/rjp5nc/vorticella
awk 'NR>1 {print $1"\t"$7}' coverage_v2/per_contig_coverage.tsv > coverage_v2/meandepth.tsv