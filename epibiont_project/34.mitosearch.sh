#!/usr/bin/env bash
#
#SBATCH -J rrna_barrnap110
#SBATCH --ntasks-per-node=16
#SBATCH -N 1
#SBATCH -t 12:00:00
#SBATCH --mem 40G
#SBATCH -o /scratch/rjp5nc/err/rrna_barrnap110.%A_%a.out
#SBATCH -e /scratch/rjp5nc/err/rrna_barrnap110.%A_%a.err
#SBATCH -p standard
#SBATCH --account berglandlab_standard
#SBATCH --mail-type=END
#SBATCH --mail-user=rjp5nc@virginia.edu


  cd /scratch/rjp5nc/vorticella
module load blast seqkit
mkdir -p vort_mito && cd vort_mito

# Ciliate mitochondrial proteins (Tetrahymena thermophila, Paramecium tetraurelia)
wget -O ciliate_mito_prot.fa \
  "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_003029,NC_001324&rettype=fasta_cds_aa&retmode=text"
grep -c ">" ciliate_mito_prot.fa

# Search all contigs >= 300 bp, ciliate genetic code (table 4 = mito code used by ciliates: UGA = Trp)
seqkit seq -m 300 ../megahit_out2/final.contigs.fa > contigs_300.fa
makeblastdb -in contigs_300.fa -dbtype nucl -out contigs_300 > /dev/null
tblastn -query ciliate_mito_prot.fa -db contigs_300 -db_gencode 4 \
  -evalue 1e-5 -num_threads 8 -max_target_seqs 20 \
  -outfmt "6 qseqid sseqid pident length evalue bitscore" > ciliate_mito_hits.tsv

# Contigs hit, with how many different mito proteins, plus depth and GC
awk -F'\t' '{g[$2][$1]=1} END {for (c in g) {n=0; for (p in g[c]) n++; print c"\t"n}}' ciliate_mito_hits.tsv \
  | sort -k2,2nr > hit_counts.tsv
awk -F'\t' 'NR==FNR {d[$1]=$2"\t"$3"\t"$4; next} {print $0"\t"d[$1]}' \
  ../coverage_v2/contig_length_depth_gc.tsv hit_counts.tsv \
  | awk 'BEGIN{print "contig\tn_proteins\tlength\tdepth\tgc"} {print}' | column -t | head -40