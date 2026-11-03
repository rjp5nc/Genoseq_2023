#!/usr/bin/env bash
#
#SBATCH -J rrna_infernal
#SBATCH --ntasks-per-node=40
#SBATCH -N 1
#SBATCH -t 48:00:00
#SBATCH --mem 60G
#SBATCH -o /scratch/rjp5nc/err/rrna_infernal.%A_%a.out
#SBATCH -e /scratch/rjp5nc/err/rrna_infernal.%A_%a.err
#SBATCH -p standard
#SBATCH --account berglandlab_standard
#SBATCH --mail-type=END
#SBATCH --mail-user=rjp5nc@virginia.edu

set -euo pipefail
module load miniforge
source activate mito

wd=/scratch/rjp5nc/vorticella
out=$wd/rrna_v2_infernal
rfam=/scratch/rjp5nc/rfamDB
silva=/scratch/rjp5nc/silvaDB/silva_ssu_nr99_dna.fa
mkdir -p "$out"
cd "$out"

if [ ! -s "$silva" ]; then
  zcat /scratch/rjp5nc/silvaDB/silva_ssu_nr99_dna.fa.gz > "$silva"
fi

seqkit seq -i "$wd"/megahit_out2/final.contigs.fa > contigs.fa
seqkit faidx contigs.fa

# Database size for E-values: total Mb x 2 (both strands)
Z=$(seqkit stats -T contigs.fa | awk 'NR==2 {printf "%.0f", $5*2/1e6}')
echo "Z = $Z Mb"

# 1. Search with Rfam rRNA models
cmscan --cpu 40 -Z "$Z" --cut_ga --rfam --nohmmonly \
  --fmt 2 --clanin "$rfam"/Rfam.clanin --oclan \
  --tblout rrna.tblout \
  "$rfam"/rrna.cm contigs.fa > rrna.cmscan.out
echo "cmscan done: $(date)"

# 2. SSU hits only, not overlapped by a better hit (olp column != "=")
#    tblout fmt2 columns: 2=model 4=contig 10=seq_from 11=seq_to 12=strand 17=score 20=olp
awk 'BEGIN{OFS="\t"} !/^#/ && $2 ~ /^SSU_/ && $20 != "=" {
       s=($10<$11)?$10:$11; e=($10<$11)?$11:$10
       print $4, s-1, e, $4"|"$2"|"s"-"e"|"$12, $17, $12 }' rrna.tblout \
  | sort -k1,1 -k2,2n -u > ssu.bed

echo "SSU hits: $(wc -l < ssu.bed)"
cut -f4 ssu.bed | cut -d'|' -f2 | sort | uniq -c

# All rRNA hits by model (incl. LSU), for reference
echo "== All rRNA models hit =="
awk '!/^#/ && $20 != "=" {print $2}' rrna.tblout | sort | uniq -c

bedtools getfasta -fi contigs.fa -bed ssu.bed -s -name > ssu_raw.fa
seqkit seq -m 300 ssu_raw.fa > ssu.fa
echo "SSU sequences >=300 bp: $(grep -c '>' ssu.fa)"

# 3. Classify against SILVA
vsearch --usearch_global ssu.fa --db "$silva" \
  --id 0.70 --threads 40 --maxaccepts 50 --top_hits_only \
  --userout silva_hits.tsv \
  --userfields query+target+id+alnlen+qcov
echo "vsearch done: $(date)"

# 4. Best hit + taxonomy + contig depth/GC/length
if [ ! -s ../rrna_v2/silva_tax.tsv ]; then
  awk '/^>/{id=substr($1,2); $1=""; sub(/^ /,""); print id"\t"$0}' "$silva" > silva_tax.tsv
else
  cp ../rrna_v2/silva_tax.tsv .
fi

awk -F'\t' 'BEGIN{OFS="\t"}
  FILENAME==ARGV[1] { tax[$1]=$2; next }
  FILENAME==ARGV[2] { if (FNR>1) { len[$1]=$2; dep[$1]=$3; gc[$1]=$4 }; next }
  { q=$1; if (q in seen) next; seen[q]=1
    split(q, p, "|"); c=p[1]
    print c, p[2], $3, $5, dep[c], gc[c], len[c], tax[$2] }' \
  silva_tax.tsv "$wd"/coverage_v2/contig_length_depth_gc.tsv silva_hits.tsv \
  | sort -t$'\t' -k5,5gr \
  | awk 'BEGIN{OFS="\t"; print "contig","model","pct_id","qcov","contig_depth","contig_gc","contig_len","silva_taxonomy"} {print}' \
  > ssu_classified.tsv

echo "Classified: $(($(wc -l < ssu_classified.tsv) - 1))"
echo "== Groups =="
awk -F'\t' 'NR>1 {split($8,t,";"); g=t[1]
  if ($8 ~ /Chloroplast/) g="Plastid"; else if ($8 ~ /Mitochondria/) g="Mitochondria"
  n[g]++} END {for (k in n) print n[k], k}' ssu_classified.tsv | sort -k1,1nr