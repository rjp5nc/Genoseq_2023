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

# Paths
HIFIASM=/scratch/rjp5nc/Reference_genomes/post_kraken/assembly.hap2_onlydaps.fasta
SPADES=/scratch/rjp5nc/spades/spades_output/contigs.fasta
HIFI=/scratch/rjp5nc/HMW/HMWDNAElvis3/m84128_250121_222443_s2.hifi_reads.bc2104.fq.gz
R1=/scratch/rjp5nc/HMW/shortreadElvis/merged_R1.fq.gz
R2=/scratch/rjp5nc/HMW/shortreadElvis/merged_R2.fq.gz

THREADS=38   # change as needed

# Filter Spades contigs >=1kb
seqkit seq -m 1000 "$SPADES" > spades.1kb.fasta

# Patch Hifiasm assembly with Spades contigs
ragtag.py patch -o ragtag_patch "$HIFIASM" spades.1kb.fasta
MINPUT=ragtag_patch/ragtag.patch.fasta

# Map HiFi reads
minimap2 -t $THREADS -x map-hifi "$MINPUT" "$HIFI" | samtools sort -@ 8 -o hifi.bam
samtools index hifi.bam

# Make fofn for NextPolish
echo "$(realpath hifi.bam)" > hifi.fofn

# Round 1: polish with HiFi
nextPolish run -g "$MINPUT" -t $THREADS -p pb -hifi_fofn hifi.fofn -o np_round1.fasta

# Map Illumina reads
minimap2 -t $THREADS -ax sr np_round1.fasta "$R1" "$R2" | samtools sort -@ 8 -o ilmn.bam
samtools index ilmn.bam

# Make fofn for Illumina
echo "$(realpath ilmn.bam)" > ilmn.fofn

# Round 2: polish with Illumina
nextPolish run -g np_round1.fasta -t $THREADS -p s -sgs_fofn ilmn.fofn -o np_round2.fasta

# Run Quast
quast -t $THREADS -o quast_np np_round2.fasta