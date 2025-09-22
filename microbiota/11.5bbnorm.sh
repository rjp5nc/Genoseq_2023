#!/usr/bin/env bash
#
#SBATCH -J bbnorm # A single job name for the array
#SBATCH --ntasks-per-node=10 # multi core
#SBATCH -N 1 # on one node
#SBATCH -t 6-00:00 # 6 days
#SBATCH --mem 40G
#SBATCH -o /scratch/rjp5nc/err/Spades%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/err/Spades%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab


#conda create -n viral_assembly -c bioconda -c conda-forge bbmap
conda activate viral_assembly

bbnorm.sh in1=/scratch/rjp5nc/UK2022_2024/unmapped_fastqs_newseq/fastqs/unmapped_trimmedmerged1.fq.gz \
         in2=/scratch/rjp5nc/UK2022_2024/unmapped_fastqs_newseq/fastqs/unmapped_trimmedmerged2.fq.gz \
         out1=/scratch/rjp5nc/UK2022_2024/unmapped_fastqs_newseq/fastqs/unmapped_norm1.fq.gz \
         out2=/scratch/rjp5nc/UK2022_2024/unmapped_fastqs_newseq/fastqs/unmapped_norm2.fq.gz \
         target=150 min=3 threads=10