#!/usr/bin/env bash
#
#SBATCH -J liftover # A single job name for the array
#SBATCH --ntasks-per-node=10 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 2-10:00 # 10 hours
#SBATCH --mem 100G
#SBATCH -o /scratch/rjp5nc/erroroutputs/lift.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/erroroutputs/lift.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab


#conda create -n plotsr
conda activate plotsr
#conda install numpy=1.21.2 pandas=1.2.4 matplotlib=3.3.4 setuptools
#conda install syri

#syri -c A_B.bam -r A.fa -q B.fa -F B --prefix A_B &

module load samtools


REF="eu_pulex_totalHiCwithallbestgapclosed.clean.fa"
QUERY="assembly.hap2_onlydaps.fasta"

cd /scratch/rjp5nc/Reference_genomes/liftover/

module load gcc/11.4.0
module load minimap2/2.28




minimap2 -x asm5 reference_genome1.fasta query_genome2.fasta > alignments.paf

paf2chain /scratch/rjp5nc/Reference_genomes/liftover/ribbon/euobtusa_vs_eupulex.paf /scratch/rjp5nc/Reference_genomes/post_kraken/${REF} /scratch/rjp5nc/Reference_genomes/post_kraken/${QUERY} /scratch/rjp5nc/Reference_genomes/chain_output.chain

chainProcessor /scratch/rjp5nc/Reference_genomes/liftover/ribbon/euobtusa_vs_eupulex.paf /scratch/rjp5nc/Reference_genomes/post_kraken/${REF} /scratch/rjp5nc/Reference_genomes/post_kraken/${QUERY} /scratch/rjp5nc/Reference_genomes/chain_output.chain

lastz /scratch/rjp5nc/Reference_genomes/post_kraken/${REF} /scratch/rjp5nc/Reference_genomes/post_kraken/${QUERY} --output=/scratch/rjp5nc/Reference_genomes/output_alignment.mfa







module load bedtools

minimap2 -ax asm5 -t 24 --eqx /scratch/rjp5nc/Reference_genomes/post_kraken/${REF} /scratch/rjp5nc/Reference_genomes/post_kraken/${QUERY} > /scratch/rjp5nc/Reference_genomes/out.sam
module load gcc/11.4.0 samtools/1.17 
samtools view -b /scratch/rjp5nc/Reference_genomes/out.sam > /scratch/rjp5nc/Reference_genomes/out.bam

samtools view -b /scratch/rjp5nc/Reference_genomes/out.bam | bedtools bamtobed > /scratch/rjp5nc/Reference_genomes/out.bed

samtools view input.bam | bam2paf > output.paf

liftOver /scratch/rjp5nc/Reference_genomes/out.bed hg19ToHg38.over.chain.gz /scratch/rjp5nc/Reference_genomes/output.bed /scratch/rjp5nc/Reference_genomes/unlifted.bed
