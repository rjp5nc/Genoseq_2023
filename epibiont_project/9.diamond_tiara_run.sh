#!/usr/bin/env bash
#
#SBATCH -J mito_search
#SBATCH --ntasks-per-node=40
#SBATCH -N 1
#SBATCH -t 12:00:00
#SBATCH --mem 60G
#SBATCH -o /scratch/rjp5nc/err/mito_search.%A_%a.out
#SBATCH -e /scratch/rjp5nc/err/mito_search.%A_%a.err
#SBATCH -p standard
#SBATCH --account berglandlab_standard
#SBATCH --mail-type=END
#SBATCH --mail-user=rjp5nc@virginia.edu

set -euo pipefail
module load miniforge
source activate mito
module load seqkit

wd=/scratch/rjp5nc/vorticella
asm=$wd/megahit_out2/final.contigs.fa
out=$wd/mito_v2
mkdir -p "$out"
cd "$out"

# 1. Tiara: composition-based classification (contigs >= 1 kb)
tiara -i "$asm" -o tiara.txt -m 1000 -t 40 --tf mit --pr
echo "Tiara done: $(date)"

# 2. DIAMOND vs RefSeq mitochondrial proteins
seqkit seq -m 1000 "$asm" > contigs_1kb.fa
diamond blastx -q contigs_1kb.fa -d /scratch/rjp5nc/mitoDB/refseq_mito \
  -o diamond_mito.tsv -p 40 --sensitive --evalue 1e-10 --max-target-seqs 1 \
  --outfmt 6 qseqid sseqid pident length evalue bitscore stitle
echo "DIAMOND done: $(date)"

# 3. Combine: one row per candidate contig
awk -F'\t' 'BEGIN{OFS="\t"}
  FILENAME==ARGV[1] { split($1, a, " "); id=a[1]; tclass[id]=$2; next }
  FILENAME==ARGV[2] { if (!($1 in best) || $6>bs[$1]) {
                        bs[$1]=$6; pid[$1]=$3
                        org=$7; sub(/.*\[/, "", org); sub(/\].*/, "", org); best[$1]=org
                        prot=$7; sub(/ \[.*/, "", prot); sub(/^[^ ]+ /, "", prot); pr[$1]=prot }
                      nhit[$1]++; next }
  FNR==1 { print $0, "tiara_class", "mito_hits", "best_pident", "best_protein", "best_organism"; next }
  { id=$1
    tc=(id in tclass)?tclass[id]:"NA"
    if (tc=="mitochondrion" || (id in best))
      print $0, tc, (id in nhit)?nhit[id]:0, (id in pid)?pid[id]:"NA",
            (id in pr)?pr[id]:"NA", (id in best)?best[id]:"NA" }' \
  tiara.txt diamond_mito.tsv "$wd"/coverage_v2/contig_length_depth_gc.tsv \
  > mito_candidates.tsv

echo "Candidates: $(($(wc -l < mito_candidates.tsv) - 1))"