#!/usr/bin/env bash
#
#SBATCH -J mitofinder_megahit
#SBATCH --ntasks-per-node=40
#SBATCH -N 1
#SBATCH -t 1:00:00
#SBATCH --mem 40G
#SBATCH -o /scratch/rjp5nc/err/mitofinder_megahit.%A_%a.out
#SBATCH -e /scratch/rjp5nc/err/mitofinder_megahit.%A_%a.err
#SBATCH -p standard
#SBATCH --account berglandlab_standard
#SBATCH --mail-type=END
#SBATCH --mail-user=rjp5nc@virginia.edu


cd /scratch/rjp5nc/vorticella/mitofinder
module load seqkit blast

# 1. Reference mitochondrial proteins (all 15 ciliate mitogenomes, incl. Thuricola)
cp test_run/translated_genes_for_database.fasta ref_mito_proteins.faa
grep -c ">" ref_mito_proteins.faa

# 2. Search the full MEGAHIT assembly (all contigs >= 300 bp), mito code 4
seqkit seq -m 300 megahit_all.fa > megahit_300.fa
makeblastdb -in megahit_300.fa -dbtype nucl -out megahit_300 > /dev/null
tblastn -query ref_mito_proteins.faa -db megahit_300 -db_gencode 4 \
  -evalue 1e-5 -num_threads 16 -max_target_seqs 50 \
  -outfmt "6 qseqid sseqid pident length evalue bitscore" > mito_tblastn.tsv
wc -l mito_tblastn.tsv