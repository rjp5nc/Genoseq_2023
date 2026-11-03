#!/usr/bin/env bash
#
#SBATCH -J hifi_check
#SBATCH --ntasks-per-node=16
#SBATCH -N 1
#SBATCH -t 12:00:00
#SBATCH --mem 100G
#SBATCH -o /scratch/rjp5nc/err/hifi_check.%A_%a.out
#SBATCH -e /scratch/rjp5nc/err/hifi_check.%A_%a.err
#SBATCH -p standard
#SBATCH --account berglandlab_standard
#SBATCH --mail-type=END
#SBATCH --mail-user=rjp5nc@virginia.edu

set -euo pipefail
module purge
module load gcc/14.2.0 minimap2/2.30 samtools seqkit 

f=/scratch/rjp5nc/HMW/HMWDNAElvis3/m84128_250121_222443_s2.hifi_reads.bc2104.fq.gz
cd /scratch/rjp5nc/vorticella/hifi_check

# Map HiFi reads; keep primary alignments with MAPQ>=10
minimap2 -t 16 -x map-hifi --secondary=no vort_ref.fa "$f" > hifi_vs_vort.paf

# PAF columns: 1 read, 2 read_len, 6 target, 7 target_len, 8 t_start, 9 t_end, 10 matches, 11 aln_len, 12 MAPQ
awk -F'\t' '$12>=10' hifi_vs_vort.paf > hifi_vs_vort.q10.paf

echo "== Reads hitting each category =="
awk -F'\t' '{cat=($6 ~ /^mito_/)?"mito":($6 ~ /^rDNA_/)?"rDNA":"nuclear"; if (!seen[$1 FS cat]++) n[cat]++}
  END {for (k in n) print k"\t"n[k]}' hifi_vs_vort.q10.paf

echo "== Mito hits (read, read_len, contig, contig_len, t_start, t_end, identity) =="
awk -F'\t' '$6 ~ /^mito_/ {printf "%s\t%d\t%s\t%d\t%d\t%d\t%.3f\n", $1, $2, $6, $7, $8, $9, $10/$11}' \
  hifi_vs_vort.q10.paf | sort -k3,3 -k2,2nr

# Save any mito-hitting reads for assembly
cut -f1 <(awk -F'\t' '$6 ~ /^mito_/' hifi_vs_vort.q10.paf) | sort -u > mito_read_ids.txt
seqkit grep -f mito_read_ids.txt "$f" -o hifi_mito_reads.fq.gz
seqkit stats hifi_mito_reads.fq.gz