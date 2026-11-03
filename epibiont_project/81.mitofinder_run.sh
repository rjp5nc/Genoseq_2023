#!/usr/bin/env bash
#
#SBATCH -J mito_tblastn
#SBATCH --ntasks-per-node=16
#SBATCH -N 1
#SBATCH -t 1:00:00
#SBATCH --mem 32G
#SBATCH -o /scratch/rjp5nc/err/mito_tblastn.%A_%a.out
#SBATCH -e /scratch/rjp5nc/err/mito_tblastn.%A_%a.err
#SBATCH -p standard
#SBATCH --account berglandlab_standard
#SBATCH --mail-type=END
#SBATCH --mail-user=rjp5nc@virginia.edu

set -euo pipefail
module purge
module load seqkit blast

wd=/scratch/rjp5nc/vorticella
cd "$wd"/mitofinder

# 1. Assembly (clean IDs, >= 300 bp) and depth table from MEGAHIT headers
seqkit seq -i "$wd"/megahit_out2/final.contigs.fa | seqkit seq -m 300 > megahit_300.fa
grep ">" "$wd"/megahit_out2/final.contigs.fa \
  | sed 's/>//; s/ flag=[^ ]*//; s/ multi=/\t/; s/ len=/\t/' > megahit_depth.tsv
seqkit fx2tab -n -l -g megahit_300.fa | awk -F'\t' '{print $1"\t"$2"\t"$3}' > megahit_300_len_gc.tsv
echo "Contigs >= 300 bp: $(wc -l < megahit_300_len_gc.tsv)"

# 2. Reference mito proteins (from the MitoFinder test run)
cp test_run/translated_genes_for_database.fasta ref_mito_proteins.faa
echo "Reference proteins: $(grep -c '>' ref_mito_proteins.faa)"

# 3. tblastn, mitochondrial code 4
makeblastdb -in megahit_300.fa -dbtype nucl -out megahit_300 > /dev/null
tblastn -query ref_mito_proteins.faa -db megahit_300 -db_gencode 4 \
  -evalue 1e-5 -num_threads 16 -max_target_seqs 50 \
  -outfmt "6 qseqid sseqid pident length evalue bitscore" > mito_tblastn.tsv
echo "tblastn hits: $(wc -l < mito_tblastn.tsv)"

# 4. Summary per contig: number of genes hit, length, depth, GC, best identity, gene names
awk -F'\t' '
  FILENAME==ARGV[1] { dep[$1]=$2; next }
  FILENAME==ARGV[2] { len[$1]=$2; gc[$1]=$3; next }
  { g=$1; sub(/^.*\|/, "", g); sub(/_[0-9]+$/, "", g)
    if (!(($2,g) in seen)) { seen[$2,g]=1; genes[$2]=genes[$2] (genes[$2]?",":"") g; n[$2]++ }
    if ($3>best[$2]) best[$2]=$3 }
  END { for (c in n) printf "%s\t%d\t%s\t%s\t%s\t%.1f\t%s\n", c, n[c], len[c], dep[c], gc[c], best[c], genes[c] }' \
  megahit_depth.tsv megahit_300_len_gc.tsv mito_tblastn.tsv \
  | sort -t$'\t' -k2,2nr \
  | awk 'BEGIN{OFS="\t"; print "contig","n_genes","length","depth","gc","best_pid","genes"} {print}' \
  > mito_contigs_summary.tsv

echo "== Vorticella-like (GC<36, depth 60-300) =="
awk -F'\t' 'NR==1 || ($5<36 && $4>=60 && $4<=300)' mito_contigs_summary.tsv
echo "== Second-ciliate-like (GC<36, depth 25-60) =="
awk -F'\t' 'NR==1 || ($5<36 && $4>=25 && $4<60)' mito_contigs_summary.tsv