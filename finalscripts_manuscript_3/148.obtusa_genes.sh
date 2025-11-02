#!/usr/bin/env bash

#SBATCH -J BEAST # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 6-0:00:00 ### 15 seconds
#SBATCH --mem 60G
#SBATCH -o /scratch/rjp5nc/erroroutputs/beast.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/erroroutputs/beast.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab
#SBATCH --mail-type=END               # Send email at job completion
#SBATCH --mail-user=rjp5nc@virginia.edu    # Email address for notifications

module load samtools





cd /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/Daphnia_obtusa_genes

/scratch/rjp5nc/Reference_genomes/post_kraken/US_obtusa_onlydaps

cut -f1 /scratch/rjp5nc/Reference_genomes/post_kraken/US_obtusa_onlydaps.fa.fai | head -n 12 > contigs12.txt

seqkit head -n 12 /scratch/rjp5nc/Reference_genomes/post_kraken/US_obtusa_onlydaps.fa > first12_contigs.fasta
samtools faidx Daphnia_obtusa_FS6_genome.fa
cut -f1 Daphnia_obtusa_FS6_genome.fa.fai | head -n 12 > contigs12.txt

cut -f1,2 Daphnia_obtusa_FS6_genome.fa.fai 
cut -f1,2 first12_contigs.fasta.fai


module load bcftools
module load samtools
module load htslib

bcftools view \
  -r JAACYE010000001.1,JAACYE010000002.1,JAACYE010000003.1,JAACYE010000004.1,JAACYE010000005.1,JAACYE010000006.1,JAACYE010000007.1,JAACYE010000008.1,JAACYE010000009.1,JAACYE010000010.1,JAACYE010000011.1,JAACYE010000012.1 \
  -Oz -o trimmed10bp_filtered_Gilmer_12.vcf.gz \
  trimmed10bp_filtered_Gilmer.vcf.gz

#zcat trimmed10bp_filtered_Gilmer_12.vcf.gz \
  | sed 's/JAACYE010000001.1/DOBTUS_1/g; 
         s/JAACYE010000002.1/DOBTUS_2/g; 
         s/JAACYE010000003.1/DOBTUS_3/g; 
         s/JAACYE010000004.1/DOBTUS_4/g; 
         s/JAACYE010000005.1/DOBTUS_5/g; 
         s/JAACYE010000006.1/DOBTUS_6/g; 
         s/JAACYE010000007.1/DOBTUS_7/g; 
         s/JAACYE010000008.1/DOBTUS_8/g; 
         s/JAACYE010000009.1/DOBTUS_9/g; 
         s/JAACYE010000010.1/DOBTUS_10/g; 
         s/JAACYE010000011.1/DOBTUS_11/g; 
         s/JAACYE010000012.1/DOBTUS_12/g' \
  | bgzip > trimmed10bp_filtered_Gilmer_12_renamed.vcf.gz

# Index
tabix -p vcf trimmed10bp_filtered_Gilmer_12_renamed.vcf.gz


cut -f1 trimmed10bp_filtered_Gilmer_12_renamed.vcf.gz | head -n 1
