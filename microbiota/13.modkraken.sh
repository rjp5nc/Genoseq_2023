#!/usr/bin/env bash
#
#SBATCH -J run_kraken # A single job name for the array
#SBATCH --ntasks-per-node=10 # multi core
#SBATCH -N 1 # on one node
#SBATCH -t 3-03:00 # 3 hours
#SBATCH --mem 100G
#SBATCH -o /scratch/rjp5nc/err/cilliates%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/err/cilliates%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

module load kraken2
export KRAKEN2_DATA_PATH="/scratch/rjp5nc/krakenDB/nt"

cd /scratch/rjp5nc/krakenDB/cilliates

esearch -db nucleotide -query "Vorticella[Organism]" | efetch -format fasta > vorticella_only.fasta

# --- Add sequences to the Kraken2 library ---
kraken2-build --add-to-library vorticella_only.fasta --db /scratch/rjp5nc/krakenDB/nt
kraken2-build --add-to-library cilliates.fasta --db /scratch/rjp5nc/krakenDB/nt
# --- Rebuild the Kraken2 database ---
kraken2-build --build --db /scratch/rjp5nc/krakenDB/nt --threads 10