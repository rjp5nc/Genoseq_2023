#!/usr/bin/env bash

#SBATCH -J crossmap # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 12:00:00 
#SBATCH --mem 20G
#SBATCH -o /scratch/rjp5nc/erroroutputs/crossmap.out # Standard output
#SBATCH -e /scratch/rjp5nc/erroroutputs/crossmap.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab
#SBATCH --mail-type=END               # Send email at job completion
#SBATCH --mail-user=rjp5nc@virginia.edu    # Email address for notifications

#Load necessary modules (if needed)
module load gcc/11.4.0
module load openmpi/4.1.4
module load python/3.11.4

#pip install --user .

#species=uspulex
#species=usobtusa
#species=usambigua
species=euobtusa
speciescrossmap=eu_obtusa
vcffile=trimmed10bp_masked_${species}.vcf.gz
sourcegenome=totalHiCwithallbestgapclosed.fa
#Daphnia_ambigua_Q001_genome.fa
#US_obtusa_onlydaps.fa
#us_pulex_ref_kap4.fa

samtools faidx /scratch/rjp5nc/Reference_genomes/post_kraken/$sourcegenome

/home/rjp5nc/.local/bin/CrossMap vcf --chromid a /scratch/rjp5nc/lastz/$speciescrossmap/chainnet/liftover.chain \
 /scratch/rjp5nc/Reference_genomes/post_kraken/$sourcegenome \
 /scratch/rjp5nc/UK2022_2024/daphnia_phylo/trimmed_10bp_repeatmasked_vcf/$vcffile \
 /scratch/rjp5nc/UK2022_2024/daphnia_phylo/trimmed_10bp_repeatmasked_vcf/lifted_${species}.vcf

#cd /home/rjp5nc/Genoseq_2023/shortread/liftover


