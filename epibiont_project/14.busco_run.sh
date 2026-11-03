#!/usr/bin/env bash
#
#SBATCH -J busco_vort
#SBATCH --ntasks-per-node=40
#SBATCH -N 1
#SBATCH -t 24:00:00
#SBATCH --mem 100G
#SBATCH -o /scratch/rjp5nc/err/busco_vort.%A_%a.out
#SBATCH -e /scratch/rjp5nc/err/busco_vort.%A_%a.err
#SBATCH -p standard
#SBATCH --account berglandlab_standard
#SBATCH --mail-type=END
#SBATCH --mail-user=rjp5nc@virginia.edu

set -euo pipefail
module purge
module load apptainer/1.4.5
module load busco/6.1.0

BUSCO="apptainer exec -B /scratch $CONTAINERDIR/busco-6.1.0.sif busco"

cd /scratch/rjp5nc/vorticella/vort_bin
db=/scratch/rjp5nc/buscoDB

for lin in alveolata_odb12 eukaryota_odb12; do
  # Standard genetic code
  $BUSCO -i vort_bin.fa -m genome -l "$lin" -c 40 \
    --offline --download_path "$db" \
    --metaeuk \
    -o busco_${lin}_std -f

  # Ciliate code (UAA/UAG = Gln, NCBI table 6)
  $BUSCO -i vort_bin.fa -m genome -l "$lin" -c 40 \
    --offline --download_path "$db" \
    --metaeuk \
    --metaeuk_parameters="--translation-table=6" \
    --metaeuk_rerun_parameters="--translation-table=6" \
    -o busco_${lin}_tt6 -f
done

echo "===== Summaries ====="
for d in busco_*; do
  echo "== $d =="
  grep -h "C:" "$d"/short_summary*.txt
done