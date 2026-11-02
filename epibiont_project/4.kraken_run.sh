#!/usr/bin/env bash
#
#SBATCH -J kraken2_v2
#SBATCH --ntasks-per-node=16
#SBATCH -N 1
#SBATCH -t 4:00:00
#SBATCH --mem 30G
#SBATCH -o /scratch/rjp5nc/err/kraken2_v2.%A_%a.out
#SBATCH -e /scratch/rjp5nc/err/kraken2_v2.%A_%a.err
#SBATCH -p standard
#SBATCH --account berglandlab_standard
#SBATCH --mail-type=END
#SBATCH --mail-user=rjp5nc@virginia.edu

set -euo pipefail
module load kraken2
module load seqkit

cd /scratch/rjp5nc/vorticella

# Random ~1% of pairs; same seed keeps pairs matched
if [ ! -s sub2_R1.fq.gz ]; then
  seqkit sample -p 0.01 -s 42 -j 8 nonDaphnia2_R1.fq.gz -o sub2_R1.fq.gz
  seqkit sample -p 0.01 -s 42 -j 8 nonDaphnia2_R2.fq.gz -o sub2_R2.fq.gz
fi
seqkit stats sub2_R1.fq.gz sub2_R2.fq.gz

kraken2 --db /scratch/rjp5nc/krakenDB/pluspf16 --threads 16 --paired \
  --confidence 0.1 \
  --report sub2_kraken_report.txt --output /dev/null \
  sub2_R1.fq.gz sub2_R2.fq.gz

head -5 sub2_kraken_report.txt
awk '$4=="P"' sub2_kraken_report.txt | sort -k1,1nr | head -20
grep -i ciliophora sub2_kraken_report.txt || echo "no Ciliophora hits"
awk '$1 > 1' sub2_kraken_report.txt > sub2_kraken_above1pct.txt