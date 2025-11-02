#!/usr/bin/env bash
#
#SBATCH -J Spades # A single job name for the array
#SBATCH --ntasks-per-node=30 # multi core
#SBATCH -N 1 # on one node
#SBATCH -t 6-00:00 # 6 days
#SBATCH --mem 150G
#SBATCH -o /scratch/rjp5nc/err/Spades%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/err/Spades%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

#conda create -n spades -c bioconda -c conda-forge spades
source $(conda info --base)/etc/profile.d/conda.sh
conda activate spades && \
spades.py \
    -1 /scratch/rjp5nc/HMW/shortreadElvis/merged_R1.fq.gz \
    -2 /scratch/rjp5nc/HMW/shortreadElvis/merged_R2.fq.gz \
    --pacbio /scratch/rjp5nc/HMW/HMWDNAElvis3/m84128_250121_222443_s2.hifi_reads.bc2104.fastq \
    -o /scratch/rjp5nc/spades/spades_output \
    -t $(nproc) \
    -m $(( $(grep MemTotal /proc/meminfo | awk '{print $2}') / 1024 / 1024 - 2 ))


awk '/^>/ {if (seqlen){print name"\t"seqlen}; name=$0; sub(/^>/,"",name); seqlen=0; next} {seqlen+=length($0)} END {print name"\t"seqlen}' /scratch/rjp5nc/spades/spades_output/scaffolds.fasta > /scratch/rjp5nc/spades/spades_output/scaffoldsize.txt

grep -c "^>" /scratch/rjp5nc/spades/spades_output/scaffolds.fasta

awk 'BEGIN{FS="\n"; RS=">"; OFS="\n"} 
NR>1 {
    seq=""
    for (i=2; i<=NF; i++) seq=seq $i
    if (length(seq) >= 14000) {
        print ">" $1
        for (i=2; i<=NF; i++) print $i
    }  
}' /scratch/rjp5nc/spades/spades_output/scaffolds.fasta > /scratch/rjp5nc/spades/spades_output/scaffolds_min14k.fasta

grep -c "^>" /scratch/rjp5nc/spades/spades_output/scaffolds_min14k.fasta
