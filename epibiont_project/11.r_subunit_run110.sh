#!/usr/bin/env bash
#
#SBATCH -J rrna_barrnap110
#SBATCH --ntasks-per-node=40
#SBATCH -N 1
#SBATCH -t 12:00:00
#SBATCH --mem 60G
#SBATCH -o /scratch/rjp5nc/err/rrna_barrnap110.%A_%a.out
#SBATCH -e /scratch/rjp5nc/err/rrna_barrnap110.%A_%a.err
#SBATCH -p standard
#SBATCH --account berglandlab_standard
#SBATCH --mail-type=END
#SBATCH --mail-user=rjp5nc@virginia.edu

set -euo pipefail
module load miniforge
source activate mito

wd=/scratch/rjp5nc/vorticella
out=$wd/rrna_v2_barrnap110
silvadir=/scratch/rjp5nc/silvaDB
silva=$silvadir/silva_ssu_nr99_dna.fa
mkdir -p "$out"
cd "$out"

barrnap --version

if [ ! -s "$silva" ]; then
  zcat "$silvadir"/silva_ssu_nr99_dna.fa.gz > "$silva"
fi

seqkit seq -i "$wd"/megahit_out2/final.contigs.fa > contigs.fa
seqkit faidx contigs.fa

# 1. Find rRNA genes (1.10.5 kingdoms: bac, arc, fun)
for k in bac arc fun; do
  barrnap --kingdom "$k" --threads 40 contigs.fa \
    > barrnap_"$k".gff 2> barrnap_"$k".log
  echo "== $k: $(grep -vc '^#' barrnap_"$k".gff) features =="
  grep -v '^#' barrnap_"$k".gff | head -3 || true
done
echo "barrnap done: $(date)"

# 2. Keep small-subunit genes only (robust to naming), exclude LSU/5S/5.8S
cat barrnap_bac.gff barrnap_arc.gff barrnap_fun.gff \
  | awk -F'\t' 'BEGIN{OFS="\t"; IGNORECASE=1}
      !/^#/ && NF>=9 &&
      $9 ~ /(16S|18S|SSU|small[ _]subunit)/ &&
      $9 !~ /(23S|28S|25S|26S|5S_|5\.8S|LSU|large[ _]subunit)/ {
        nm=$9
        if (match(nm, /Name=[^;]+/)) nm=substr(nm, RSTART+5, RLENGTH-5)
        else if (match(nm, /product=[^;]+/)) nm=substr(nm, RSTART+8, RLENGTH-8)
        gsub(/[ |]/, "_", nm)
        print $1, $4-1, $5, $1"|"nm"|"$4"-"$5"|"$7, $6, $7 }' \
  | sort -k1,1 -k2,2n -u > ssu.bed

echo "SSU genes found: $(wc -l < ssu.bed)"
cut -f4 ssu.bed | cut -d'|' -f2 | sort | uniq -c

bedtools getfasta -fi contigs.fa -bed ssu.bed -s -name > ssu_raw.fa
seqkit seq -m 300 ssu_raw.fa > ssu.fa
echo "SSU sequences >=300 bp: $(grep -c '>' ssu.fa || true)"

# 3. Classify against SILVA
vsearch --usearch_global ssu.fa --db "$silva" \
  --id 0.70 --threads 40 --maxaccepts 50 --top_hits_only \
  --userout silva_hits.tsv \
  --userfields query+target+id+alnlen+qcov
echo "vsearch done: $(date)"

# 4. Best hit + taxonomy + contig depth/GC/length
if [ -s ../rrna_v2_infernal/silva_tax.tsv ]; then
  cp ../rrna_v2_infernal/silva_tax.tsv .
else
  awk '/^>/{id=substr($1,2); $1=""; sub(/^ /,""); print id"\t"$0}' "$silva" > silva_tax.tsv
fi

awk -F'\t' 'BEGIN{OFS="\t"}
  FILENAME==ARGV[1] { tax[$1]=$2; next }
  FILENAME==ARGV[2] { if (FNR>1) { len[$1]=$2; dep[$1]=$3; gc[$1]=$4 }; next }
  { q=$1; if (q in seen) next; seen[q]=1
    split(q, p, "|"); c=p[1]
    print c, p[2], $3, $5, dep[c], gc[c], len[c], tax[$2] }' \
  silva_tax.tsv "$wd"/coverage_v2/contig_length_depth_gc.tsv silva_hits.tsv \
  | sort -t$'\t' -k5,5gr \
  | awk 'BEGIN{OFS="\t"; print "contig","gene","pct_id","qcov","contig_depth","contig_gc","contig_len","silva_taxonomy"} {print}' \
  > ssu_classified.tsv

echo "Classified: $(($(wc -l < ssu_classified.tsv) - 1))"
echo "== Groups =="
awk -F'\t' 'NR>1 {split($8,t,";"); g=t[1]
  if ($8 ~ /Chloroplast/) g="Plastid"; else if ($8 ~ /Mitochondria/) g="Mitochondria"
  n[g]++} END {for (k in n) print n[k], k}' ssu_classified.tsv | sort -k1,1nr