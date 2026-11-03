#!/usr/bin/env bash
#
#SBATCH -J kraken_big
#SBATCH --ntasks-per-node=16
#SBATCH -N 1
#SBATCH -t 1:00:00
#SBATCH --mem 40G
#SBATCH -o /scratch/rjp5nc/err/kraken_big.%A_%a.out
#SBATCH -e /scratch/rjp5nc/err/kraken_big.%A_%a.err
#SBATCH -p standard
#SBATCH --account berglandlab_standard


set -euo pipefail
module purge
module load kraken2 seqkit

cd /scratch/rjp5nc/vorticella
asm=megahit_out2/final.contigs.fa

# Contigs > 70 kb (rebuild if missing)
if [ ! -s contigs_over70kb.fa ] || [ ! -s contigs_over70kb.tsv ]; then
  grep ">" $asm | sed 's/>//' \
    | awk '{m=$3; l=$4; sub("multi=","",m); sub("len=","",l); if (l+0>70000) print $1"\t"m"\t"l}' \
    | sort -k3,3nr > contigs_over70kb.tsv
  seqkit grep -f <(cut -f1 contigs_over70kb.tsv) $asm | seqkit seq -i > contigs_over70kb.fa
fi
seqkit fx2tab -n -i -g contigs_over70kb.fa > contigs_over70kb_gc.tsv
echo "Contigs > 70 kb: $(grep -c '>' contigs_over70kb.fa)"

# Copy DB to node-local disk for fast loading
localdb=${TMPDIR:-/tmp}/pluspf16
mkdir -p "$localdb"
cp /scratch/rjp5nc/krakenDB/pluspf16/*.k2d "$localdb"/

# Kraken without a confidence threshold
kraken2 --db "$localdb" --threads 16 --use-names \
  --report contigs_over70kb_kraken_noconf_report.txt \
  --output contigs_over70kb_kraken_noconf.out \
  contigs_over70kb.fa
rm -rf "$localdb"

# Table: contig, length, coverage, GC, Kraken assignment
awk -F'\t' '
  FILENAME==ARGV[1] { cov[$1]=$2; len[$1]=$3; next }
  FILENAME==ARGV[2] { gc[$1]=$2; next }
  { print $2"\t"len[$2]"\t"cov[$2]"\t"gc[$2]"\t"$3 }' \
  contigs_over70kb.tsv contigs_over70kb_gc.tsv contigs_over70kb_kraken_noconf.out \
  | sort -t$'\t' -k2,2nr \
  | awk 'BEGIN{OFS="\t"; print "contig","length","coverage","gc","kraken"} {print}' \
  > contigs_over70kb_classified_noconf.tsv

column -t -s$'\t' contigs_over70kb_classified_noconf.tsv

echo "== Summary by taxon =="
awk -F'\t' 'NR>1 {n[$5]++; s[$5]+=$2} END {for (k in n) printf "%d contigs\t%.2f Mb\t%s\n", n[k], s[k]/1e6, k}' \
  contigs_over70kb_classified_noconf.tsv | sort -t$'\t' -k2,2nr

echo "== Higher-level breakdown (from report, >=1% of contigs) =="
awk -F'\t' '$1+0>=1 && ($4=="D" || $4=="P" || $4=="C" || $4=="G")' contigs_over70kb_kraken_noconf_report.txt