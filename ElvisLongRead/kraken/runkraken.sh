#!/usr/bin/env bash
#
#SBATCH -J run_kraken # A single job name for the array
#SBATCH --ntasks-per-node=30 # multi core
#SBATCH -N 1 # on one node
#SBATCH -t 1-03:00 # 3 hours
#SBATCH --mem 150G
#SBATCH -o /scratch/rjp5nc/Canu_error/kraken%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/Canu_error/kraken%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab


module load kraken2

KRAKEN2_DATA_PATH="/scratch/rjp5nc/krakenDB/nt"

echo 'export KRAKEN2_DATA_PATH="/scratch/rjp5nc/krakenDB/nt"' >> ~/.bashrc
source ~/.bashrc

kraken2 --memory-mapping --db /scratch/rjp5nc/krakenDB/nt --threads 30 --report /scratch/rjp5nc/krakenDB/US_pulex/report.txt --classified-out /scratch/rjp5nc/krakenDB/US_pulex/kraken_classified_output.txt --output /scratch/rjp5nc/krakenDB/US_pulex/kraken_output.txt --use-names /scratch/rjp5nc/Reference_genomes/orig_ref/GCA_900092285.2_PA42_4.1_genomic.fna

grep '^>' /scratch/rjp5nc/krakenDB/US_pulex/kraken_classified_output.txt > /scratch/rjp5nc/krakenDB/US_pulex/classified_headers.txt


#kraken2-build --standard --db /scratch/rjp5nc/krakenDB/test2