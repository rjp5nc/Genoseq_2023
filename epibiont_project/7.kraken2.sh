#!/usr/bin/env bash
#
#SBATCH -J kraken_allcontigs
#SBATCH --ntasks-per-node=16
#SBATCH -N 1
#SBATCH -t 4:00:00
#SBATCH --mem 40G
#SBATCH -o /scratch/rjp5nc/err/kraken_allcontigs.%A_%a.out
#SBATCH -e /scratch/rjp5nc/err/kraken_allcontigs.%A_%a.err
#SBATCH -p standard
#SBATCH --account berglandlab_standard
#SBATCH --mail-type=END
#SBATCH --mail-user=rjp5nc@virginia.edu

set -euo pipefail
module load kraken2

wd=/scratch/rjp5nc/vorticella
out=$wd/kraken_contigs_v2
mkdir -p "$out"

# Copy DB to node-local disk for fast loading
localdb=${TMPDIR:-/tmp}/pluspf16
mkdir -p "$localdb"
cp /scratch/rjp5nc/krakenDB/pluspf16/*.k2d "$localdb"/
echo "DB copied: $(date)"

cd "$wd"
kraken2 --db "$localdb" --threads 16 \
  --confidence 0.1 --use-names \
  --report "$out"/report.txt \
  --output "$out"/per_contig.out \
  megahit_out2/final.contigs.fa
echo "Kraken done: $(date)"
rm -rf "$localdb"

# taxid -> domain, phylum (from the report's tree structure)
awk -F'\t' '{
    name=$6; match(name, /^ */); d=RLENGTH/2; sub(/^ +/, "", name)
    rk[d]=$4; nm[d]=name
    dom="unclassified"; phy="NA"
    for (i=0; i<=d; i++) { if (rk[i]=="D") dom=nm[i]; if (rk[i]=="P") phy=nm[i] }
    if ($5==1 || $5==131567) dom="root_only"
    print $5"\t"dom"\t"phy
  }' "$out"/report.txt > "$out"/taxid_lineage.tsv

# Per-contig table: contig, length, depth, gc, kraken taxon, domain, phylum
awk -F'\t' 'BEGIN{OFS="\t"}
  FILENAME==ARGV[1] { dom[$1]=$2; phy[$1]=$3; next }
  FILENAME==ARGV[2] { t=$3; tid=0
                      if (match(t, /taxid [0-9]+/)) tid=substr(t, RSTART+6, RLENGTH-6)
                      sub(/ \(taxid [0-9]+\)/, "", t)
                      taxon[$2]=t; tax[$2]=tid; next }
  FNR==1 { print $0, "kraken_taxon", "domain", "phylum"; next }
  { tid=($1 in tax) ? tax[$1] : 0
    d=(tid==0) ? "unclassified" : ((tid in dom) ? dom[tid] : "NA")
    p=(tid==0) ? "NA" : ((tid in phy) ? phy[tid] : "NA")
    print $0, taxon[$1], d, p }' \
  "$out"/taxid_lineage.tsv "$out"/per_contig.out coverage_v2/contig_length_depth_gc.tsv \
  > "$out"/contig_table.tsv

head -5 "$out"/contig_table.tsv