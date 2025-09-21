#!/usr/bin/env bash
#
#SBATCH -J Spades # A single job name for the array
#SBATCH --ntasks-per-node=30 # multi core
#SBATCH -N 1 # on one node
#SBATCH -t 6-00:00 # 6 days
#SBATCH --mem 300G
#SBATCH -o /scratch/rjp5nc/err/Spades%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/err/Spades%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

#conda create -n spades -c bioconda -c conda-forge spades
source ~/miniconda3/etc/profile.d/conda.sh
conda activate spades && \
#spades.py --meta \
#    -1 /scratch/rjp5nc/UK2022_2024/unmapped_fastqs_newseq/fastqs/unmapped_trimmedmerged1.fq.gz \
#    -2 /scratch/rjp5nc/UK2022_2024/unmapped_fastqs_newseq/fastqs/unmapped_trimmedmerged2.fq.gz \
#    -o /scratch/rjp5nc/UK2022_2024/unmapped_fastqs_newseq/SPADES \
#    -t $(nproc) \
#    -m 240 \
#    --tmp-dir /scratch/rjp5nc/UK2022_2024/unmapped_fastqs_newseq/tmp

#-t $(nproc) 

    spades.py \
    --restart-from last \
    -k 21,33 \ 
    -o /scratch/rjp5nc/UK2022_2024/unmapped_fastqs_newseq/SPADES \
    -t $SLURM_CPUS_PER_TASK\
    -m 240 \
    --tmp-dir /scratch/rjp5nc/UK2022_2024/unmapped_fastqs_newseq/tmp

#spades.py --meta \
#  --only-assembler \
#  -1 /scratch/rjp5nc/UK2022_2024/unmapped_fastqs_newseq/SPADES/corrected/unmapped_trimmedmerged1.00.0_0.cor.fastq.gz \
##  -2 /scratch/rjp5nc/UK2022_2024/unmapped_fastqs_newseq/SPADES/corrected/unmapped_trimmedmerged2.00.0_0.cor.fastq.gz \
 # --trusted-contigs /scratch/rjp5nc/UK2022_2024/unmapped_fastqs_newseq/SPADES/K21/contigs.fasta \
 # --trusted-contigs /scratch/rjp5nc/UK2022_2024/unmapped_fastqs_newseq/SPADES/K33/contigs.fasta \
 # -o /scratch/rjp5nc/UK2022_2024/unmapped_fastqs_newseq/SPADES_retry \
 # -t $(nproc) \
 # -m 240 \
 # --tmp-dir /scratch/rjp5nc/UK2022_2024/unmapped_fastqs_newseq/tmp

#awk '/^>/ {if (seqlen){print name"\t"seqlen}; name=$0; sub(/^>/,"",name); seqlen=0; next} {seqlen+=length($0)} END {print name"\t"seqlen}' /scratch/rjp5nc/UK2022_2024/unmapped_fastqs_newseq/SPADES/scaffolds.fasta > /scratch/rjp5nc/UK2022_2024/unmapped_fastqs_newseq/SPADES/scaffoldsize.txt

#grep -c "^>" /scratch/rjp5nc/UK2022_2024/unmapped_fastqs_newseq/SPADES/scaffolds.fasta

#awk 'BEGIN{FS="\n"; RS=">"; OFS="\n"} 
#NR>1 {
#    seq=""
#    for (i=2; i<=NF; i++) seq=seq $i
 #   if (length(seq) >= 14000) {
 #       print ">" $1
 #       for (i=2; i<=NF; i++) print $i
 #   }  
#}' /scratch/rjp5nc/UK2022_2024/unmapped_fastqs_newseq/SPADES/scaffolds.fasta > /scratch/rjp5nc/UK2022_2024/unmapped_fastqs_newseq/SPADES/scaffolds_min14k.fasta

#grep -c "^>" /scratch/rjp5nc/UK2022_2024/unmapped_fastqs_newseq/SPADES/scaffolds_min14k.fasta
