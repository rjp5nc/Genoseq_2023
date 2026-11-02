#!/usr/bin/env bash
#
#SBATCH -J move # A single job name for the array
#SBATCH --ntasks-per-node=40 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 48:00:00
#SBATCH --mem 60G
#SBATCH -o /scratch/rjp5nc/err/move.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/err/move.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab_standard
#SBATCH --mail-type=END               # Send email at job completion
#SBATCH --mail-user=rjp5nc@virginia.edu    # Email address for notifications


set -euo pipefail
module load kraken2

cd /scratch/rjp5nc/vorticella

if [ ! -s sub_R1.fq ]; then
  zcat nonDaphnia_R1.fq.gz | head -4000000 > sub_R1.fq
  zcat nonDaphnia_R2.fq.gz | head -4000000 > sub_R2.fq
fi

kraken2 --db /scratch/rjp5nc/krakenDB/pluspf16 --threads 16 --paired \
  --report sub_kraken_report.txt --output /dev/null \
  sub_R1.fq sub_R2.fq


  head -5 sub_kraken_report.txt                    # unclassified vs. root
awk '$4=="P"' sub_kraken_report.txt | sort -k1,1nr | head -20   # top phyla
grep -i ciliophora sub_kraken_report.

/scratch/rjp5nc/err/move.22275796_4294967294.out