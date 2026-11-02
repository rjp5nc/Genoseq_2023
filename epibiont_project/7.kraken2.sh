
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

  cd /scratch/rjp5nc/vorticella
module load kraken2

kraken2 --db /scratch/rjp5nc/krakenDB/pluspf16 --threads 16 \
  --confidence 0.1 --use-names \
  --report vort_candidates_kraken_report.txt \
  --output vort_candidates_kraken.out \
  vort_candidates_v1.fa

# Summary: what fraction is classified, and as what
head -3 vort_candidates_kraken_report.txt
awk '$4=="D" || $4=="P"' vort_candidates_kraken_report.txt | sort -k1,1nr | head -15