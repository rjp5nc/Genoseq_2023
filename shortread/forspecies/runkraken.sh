#!/usr/bin/env bash
#
#SBATCH -J run_kraken # A single job name for the array
#SBATCH --ntasks-per-node=30 # multi core
#SBATCH -N 1 # on one node
#SBATCH -t 4-03:00 # 3 hours
#SBATCH --mem 150G
#SBATCH -o /scratch/rjp5nc/erroroutputs/kraken%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/erroroutputs/kraken%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

module load kraken2

export KRAKEN2_DATA_PATH="/scratch/rjp5nc/krakenDB/nt"
#KRAKEN2_DATA_PATH="/scratch/rjp5nc/krakenDB/nt"

#echo 'export KRAKEN2_DATA_PATH="/scratch/rjp5nc/krakenDB/nt"' >> ~/.bashrc
#source ~/.bashrc

#gunzip -c /scratch/rjp5nc/UK2022_2024/allshortreads/01.RawData/newseq/Rockpool4_F9/Rockpool4_F9_trimmedmerged1.fq.gz | head -n 10000 | gzip > /scratch/rjp5nc/UK2022_2024/allshortreads/01.RawData/newseq/Rockpool4_F9/Rockpool4_F9_trimmedmerged1_10000.fq.gz
#gunzip -c /scratch/rjp5nc/UK2022_2024/allshortreads/01.RawData/newseq/Rockpool4_F9/Rockpool4_F9_trimmedmerged2.fq.gz | head -n 10000 | gzip > /scratch/rjp5nc/UK2022_2024/allshortreads/01.RawData/newseq/Rockpool4_F9/Rockpool4_F9_trimmedmerged2_10000.fq.gz


kraken2 --memory-mapping --db /scratch/rjp5nc/krakenDB/nt \
--threads 30 --report /scratch/rjp5nc/krakenDB/shortread/report.txt \
--classified-out /scratch/rjp5nc/krakenDB/shortread/kraken_classified_output.txt \
--output /scratch/rjp5nc/krakenDB/shortread/kraken_output.txt \
--use-names \
/scratch/rjp5nc/UK2022_2024/allshortreads/01.RawData/newseq/Rockpool4_C4/Rockpool4_C4_trimmedmerged1.fq.gz \
/scratch/rjp5nc/UK2022_2024/allshortreads/01.RawData/newseq/Rockpool4_C4/Rockpool4_C4_trimmedmerged2.fq.gz

grep '^>' /scratch/rjp5nc/krakenDB/shortread/kraken_classified_output_Rockpool4_F9.txt > /scratch/rjp5nc/krakenDB/shortread/classified_headers_Rockpool4_F9.txt


#kraken2-build --standard --db /scratch/rjp5nc/krakenDB/test2