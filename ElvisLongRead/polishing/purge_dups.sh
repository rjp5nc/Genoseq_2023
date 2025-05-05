#!/bin/bash


#DIDNT USE THIS




wd="/scratch/rjp5nc/HMW/HMWDNAElvis3/hifiasm_out"
cd ${wd}

# Set input assembly file and number of threads
HIFI_ASSEMBLY="assembly.bp.hap1.p_ctg.gfa"
THREADS=20
REFERENCE_GENOME="reference_genome.fa"  # Optional, only needed for RaGOO

echo "Step 1: Run Purge_dups to remove duplications"
gfatools gfa2fa $HIFI_ASSEMBLY > assembly.p_ctg.fa
minimap2 -x asm5 -t 20 assembly.p_ctg.fa assembly.hap2.p_ctg.fa > self_align.paf
quickmerge -d delta_file.out -q assembly.p_ctg.fa -r assembly.hap2.p_ctg.fa -p self_align.paf -hco 5.0 -c 1.5 -l 500 -ml 0.0 -o merged_scaffolds.fa





awk '/^>/ {if (seqlen > 100000) print seqname "\n" seq; seqname=$0; seq=""; next} {seq = seq $0} END {if (seqlen > 100000) print seqname "\n" seq}' assembly.p_ctg.fa > large_contigs.fa
awk '/^>/ {if (seqlen > 200000) print seqname "\n" seq; seqname=$0; seq=""; next} {seq = seq $0} END {if (seqlen > 200000) print seqname "\n" seq}' assembly.p_ctg.fa > large_contigs2.fa
awk '/^>/ {if (seqname && seqlen > 500000) print seqname "\n" seq; seqname=$0; seq=""; seqlen=0; next} {seq = seq $0; seqlen += length($0)} END {if (seqlen > 500000) print seqname "\n" seq}' assembly.p_ctg.fa > large_contigs3.fa



grep -c ">" assembly.p_ctg.fa
grep -c ">" large_contigs.fa
grep -c ">" large_contigs2.fa
grep -c ">" large_contigs3.fa
grep -c ">" merged_self_align.paf.fasta
grep -c ">" assembly.hap2.p_ctg.fa
