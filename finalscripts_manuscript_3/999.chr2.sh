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
# index the fasta (creates first12_contigs.fasta.fai if not already there)
samtools faidx /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/first12_contigs.fasta

# extract the contig
samtools faidx /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/first12_contigs.fasta JAACYE010000002.1 \
  > /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/JAACYE010000002.1.fasta

mkdir /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/chr2

mv /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/JAACYE010000002.1.fasta /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/chr2


seqkit sliding --window 150 --step 150 \
  /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/chr2/JAACYE010000002.1.fasta \
  > /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/chr2/JAACYE010000002.1_split_150.fasta

  grep "^>" /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/chr2/JAACYE010000002.1_split_10kb.fasta \
  | cut -d' ' -f1 \
  | sed 's/^>//' \
  > /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/chr2/contig_names.txt
