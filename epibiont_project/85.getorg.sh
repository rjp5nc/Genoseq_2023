#!/usr/bin/env bash
#
#SBATCH -J getorganelle
#SBATCH --ntasks-per-node=16
#SBATCH -N 1
#SBATCH -t 24:00:00
#SBATCH --mem 80G
#SBATCH -o /scratch/rjp5nc/err/getorganelle.%A_%a.out
#SBATCH -e /scratch/rjp5nc/err/getorganelle.%A_%a.err
#SBATCH -p standard
#SBATCH --account berglandlab_standard
#SBATCH --mail-type=END
#SBATCH --mail-user=rjp5nc@virginia.edu

set -euo pipefail
module purge
module load seqkit
export PATH=/home/rjp5nc/miniconda3/envs/getorg/bin:$PATH
get_organelle_from_reads.py --version

wd=/scratch/rjp5nc/vorticella
out=$wd/getorganelle
mkdir -p "$out"; cd "$out"

# Seeds and labels: all low-GC mito pieces (both candidate organisms)
seqkit grep -p k141_466 -p k141_146793 -p k141_293174 -p k141_204962 -p k141_159989 \
            -p k141_359459 -p k141_346990 -p k141_129714 -p k141_277908 -p k141_151643 -p k141_335728 \
            $wd/megahit_out2/final.contigs.fa | seqkit seq -i > seeds.fa
cp seeds.fa labels.fa
echo "Seeds: $(grep -c '>' seeds.fa)"

get_organelle_from_reads.py \
  -1 $wd/nonDaphnia2_R1.fq.gz -2 $wd/nonDaphnia2_R2.fq.gz \
  -F anonym -s seeds.fa --genes labels.fa \
  -o run1 -R 15 -k 21,33,45,55,65,77 -t 16 \
  --reduce-reads-for-coverage inf --max-reads inf \
  --expected-max-size 80000 --overwrite

echo "===== Result ====="
ls run1/
seqkit fx2tab -n -l -g run1/*.path_sequence.fasta 2>/dev/null | sort -k2,2nr | head -20 || true