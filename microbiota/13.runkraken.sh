#!/usr/bin/env bash
#
#SBATCH -J run_kraken # A single job name for the array
#SBATCH --ntasks-per-node=10 # multi core
#SBATCH -N 1 # on one node
#SBATCH -t 1-03:00 # 3 hours
#SBATCH --mem 300G
#SBATCH -o /scratch/rjp5nc/err/kraken%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/err/kraken%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

module load kraken2

export KRAKEN2_DATA_PATH="/scratch/rjp5nc/krakenDB/nt"
export kraken2_DATA_PATH="/scratch/rjp5nc/krakenDB/nt"
KRAKEN2_DATA_PATH="/scratch/rjp5nc/krakenDB/nt"

echo 'export KRAKEN2_DATA_PATH="/scratch/rjp5nc/krakenDB/nt"' >> ~/.bashrc
source ~/.bashrc

#cd /scratch/rjp5nc/krakenDB/nt
#tar -xvzf k2_core_nt_20241228.tar.gz
#kraken2-build --build --threads 10 --db /scratch/rjp5nc/krakenDB/nt

cd /scratch/rjp5nc/UK2022_2024/unmapped_fastqs_newseq/SPADES_norm/

kraken2 --memory-mapping --db /scratch/rjp5nc/krakenDB/nt \
  --threads 10 \
  --use-names \
  --report scaffolds_report.txt \
  --output scaffolds_output.txt \
  /scratch/rjp5nc/UK2022_2024/unmapped_fastqs_newseq/SPADES_norm/scaffolds_min10k.fasta

#kraken2 --memory-mapping --db /scratch/rjp5nc/krakenDB/nt --threads 10 --report /scratch/rjp5nc/UK2022_2024/unmapped_fastqs_newseq/SPADES_norm/scaffolds_min10k_kraken_report.txt --classified-out /scratch/rjp5nc/UK2022_2024/unmapped_fastqs_newseq/SPADES_norm/scaffolds_min10k_kraken_classified_output.txt --output /scratch/rjp5nc/UK2022_2024/unmapped_fastqs_newseq/SPADES_norm/scaffolds_min10k_kraken_output.txt --use-names /scratch/rjp5nc/UK2022_2024/unmapped_fastqs_newseq/SPADES_norm/scaffolds_min10k.fasta

#grep '^>' /scratch/rjp5nc/UK2022_2024/unmapped_fastqs_newseq/SPADES_norm/scaffolds_min10k_kraken_classified_output.txt > /scratch/rjp5nc/UK2022_2024/unmapped_fastqs_newseq/SPADES_norm/scaffolds_min10k_classified_headers.txt


#kraken2-build --standard --db /scratch/rjp5nc/krakenDB/test2