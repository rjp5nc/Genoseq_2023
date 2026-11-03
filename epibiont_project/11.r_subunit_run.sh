#!/usr/bin/env bash
#
#SBATCH -J rrna_classify
#SBATCH --ntasks-per-node=40
#SBATCH -N 1
#SBATCH -t 12:00:00
#SBATCH --mem 60G
#SBATCH -o /scratch/rjp5nc/err/rrna_classify.%A_%a.out
#SBATCH -e /scratch/rjp5nc/err/rrna_classify.%A_%a.err
#SBATCH -p standard
#SBATCH --account berglandlab_standard
#SBATCH --mail-type=END
#SBATCH --mail-user=rjp5nc@virginia.edu

set -euo pipefail
module load miniforge
source activate mito

wd=/scratch/rjp5nc/vorticella
out=$wd/rrna_v2
silvadir=/scratch/rjp5nc/silvaDB
silva=$silvadir/silva_ssu_nr99_dna.fa
mkdir -p "$out"
cd "$out"

# SILVA: decompress once (vsearch reads plain FASTA fastest)
if [ ! -s "$silva" ]; then
  zcat "$silvadir"/silva_ssu_nr99_dna.fa.gz > "$silva"
fi

# Clean headers (contig ID only) so all tools agree on names
seqkit seq -i "$wd"/megahit_out2/final.contigs.fa > contigs.fa
seqkit faidx contigs.fa

# 1. Find rRNA genes with each kingdom model
for k in bac arc euk mito; do
  barrnap --kingdom "$k" --threads 40 --reject 0.3 contigs.fa \
    > barrnap_"$k".gff 2> barrnap_"$k".log
done
echo "barrnap done: $(date)"

# 2. Keep SSU only, per kingdom model:
#    bac/arc -> 16S | euk -> 18S | mito -> 12S (mito 16S is the LSU, so excluded)
{
  awk -F'\t' '!/^#/ && $9 ~ /Name=16S_rRNA/' barrnap_bac.gff
  awk -F'\t' '!/^#/ && $9 ~ /Name=16S_rRNA/' barrnap_arc.gff
  awk -F'\t' '!/^#/ && $9 ~ /Name=18S_rRNA/' barrnap_euk.gff
  awk -F'\t' '!/^#/ && $9 ~ /Name=12S_rRNA/' barrnap_mito.gff
} | awk -F'\t' 'BEGIN{OFS="\t"} {
      split($9, a, ";"); nm=a[1]; sub("Name=", "", nm)
      print $1, $4-1, $5, $1"|"nm"|"$4"-"$5"|"$7, $6, $7 }' \
  | sort -k1,1 -k2,2n -u > ssu.bed

echo "SSU genes found: $(wc -l < ssu.bed)"
cut -f4 ssu.bed | cut -d'|' -f2 | sort | uniq -c

bedtools getfasta -fi contigs.fa -bed ssu.bed -s -name > ssu_raw.fa
seqkit seq -m 300 ssu_raw.fa > ssu.fa
echo "SSU sequences >=300 bp: $(grep -c '>' ssu.fa)"

# 3. Classify against SILVA
vsearch --usearch_global ssu.fa --db "$silva" \
  --id 0.70 --threads 40 --maxaccepts 50 --top_hits_only \
  --userout silva_hits.tsv \
  --userfields query+target+id+alnlen+qcov
echo "vsearch done: $(date)"

# 4. Best hit per SSU + SILVA taxonomy + contig depth/GC/length
awk '/^>/{id=substr($1,2); $1=""; sub(/^ /,""); print id"\t"$0}' "$silva" > silva_tax.tsv

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

# 5. Summary by top-level group
echo "== Groups =="
awk -F'\t' 'NR>1 {split($8,t,";"); g=t[1]
  if ($8 ~ /Chloroplast/) g="Plastid"; else if ($8 ~ /Mitochondria/) g="Mitochondria"
  n[g]++} END {for (k in n) print n[k], k}' ssu_classified.tsv | sort -k1,1nr