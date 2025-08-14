#!/usr/bin/env bash
#
#SBATCH -J mergehifiasmspades # A single job name for the array
#SBATCH --ntasks-per-node=40 # multi core
#SBATCH -N 1 # on one node
#SBATCH -t 2-00:00 # 6 days
#SBATCH --mem 100G
#SBATCH -o /scratch/rjp5nc/Canu_error/hifiasm%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/Canu_error/hifiasm%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

#conda create -n hybridmerge -c bioconda -c conda-forge ragtag minimap2 samtools nextpolish busco quast seqkit

source $(conda info --base)/etc/profile.d/conda.sh
conda activate hybridmerge

cd /scratch/rjp5nc/HMW/mergeSpadesHifiasm/

HIFIASM=/scratch/rjp5nc/Reference_genomes/post_kraken/assembly.hap2_onlydaps.fasta
SPADES=/scratch/rjp5nc/spades/spades_output/contigs.fasta
HIFI=/scratch/rjp5nc/HMW/HMWDNAElvis3/m84128_250121_222443_s2.hifi_reads.bc2104.fq.gz
R1=/scratch/rjp5nc/HMW/shortreadElvis/merged_R1.fq.gz
R2=/scratch/rjp5nc/HMW/shortreadElvis/merged_R2.fq.gz

seqkit seq -m 1000 "$SPADES" > spades.1kb.fasta

ragtag.py patch -o ragtag_patch "$HIFIASM" spades.1kb.fasta
MINPUT=ragtag_patch/ragtag.patch.fasta

minimap2 -t $SLURM_CPUS_PER_TASK -x map-hifi "$MINPUT" "$HIFI" | samtools sort -@ 8 -o hifi.bam
samtools index hifi.bam
nextPolish run -g "$MINPUT" -t $SLURM_CPUS_PER_TASK -p pb -s hifi.bam -o np_round1.fasta

minimap2 -t $SLURM_CPUS_PER_TASK -ax sr np_round1.fasta "$R1" "$R2" | samtools sort -@ 8 -o ilmn.bam
samtools index ilmn.bam
nextPolish run -g np_round1.fasta -t $SLURM_CPUS_PER_TASK -p s -s ilmn.bam -o np_round2.fasta

quast -t $SLURM_CPUS_PER_TASK -o quast_np np_round2.fasta
