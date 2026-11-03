#!/usr/bin/env bash
#
#SBATCH -J vort_mito_kmers
#SBATCH --ntasks-per-node=8
#SBATCH -N 1
#SBATCH -t 2:00:00
#SBATCH --mem 120G
#SBATCH --array=1-4
#SBATCH -o /scratch/rjp5nc/err/vort_mito_kmers.%A_%a.out
#SBATCH -e /scratch/rjp5nc/err/vort_mito_kmers.%A_%a.err
#SBATCH -p standard
#SBATCH --account berglandlab_standard
#SBATCH --mail-type=END
#SBATCH --mail-user=rjp5nc@virginia.edu

set -euo pipefail
base=/scratch/rjp5nc/vorticella/vort_mito_asm
cd "$base"

# Task 1 -> k=25, 2 -> 39, 3 -> 45, 4 -> 55
kmers=(25 39 45 55)
k=${kmers[$((SLURM_ARRAY_TASK_ID - 1))]}

run=Vort_mito_k${k}
outdir=$base/runs/$run
mkdir -p "$outdir"

novo=$(ls NOVOPlasty/NOVOPlasty*.pl | head -1)
echo "Task $SLURM_ARRAY_TASK_ID: k=$k, using $novo"

cat > "$outdir"/config.txt << EOF
Project:
-----------------------
Project name          = $run
Type                  = mito
Genome Range          = 25000-70000
K-mer                 = $k
Max memory            = 100
Extended log          = 1
Save assembled reads  = no
Seed Input            = $base/seed.fa
Extend seed directly  = no
Reference sequence    = 
Variance detection    = 
Chloroplast sequence  = 

Dataset 1:
-----------------------
Read Length           = 150
Insert size           = 300
Platform              = illumina
Single/Paired         = PE
Combined reads        = 
Forward reads         = /scratch/rjp5nc/vorticella/nonDaphnia2_R1.fq.gz
Reverse reads         = /scratch/rjp5nc/vorticella/nonDaphnia2_R2.fq.gz
Store Hash            = 

Heteroplasmy:
-----------------------
MAF                   = 
HP exclude list       = 
PCR-free              = 

Optional:
-----------------------
Insert size auto      = yes
Use Quality Scores    = no
Reduce ambigious N's  = 
Output path           = $outdir/
EOF

echo "Start: $(date)"
perl "$novo" -c "$outdir"/config.txt
echo "End: $(date)"

grep -h -A4 "Assembly 1 finished" "$outdir"/log_${run}.txt || echo "no assembly summary for k=$k"




# cd /scratch/rjp5nc/vorticella
# mkdir -p vort_mito_asm && cd vort_mito_asm

# # Get NOVOPlasty
# git clone https://github.com/ndierckx/NOVOPlasty.git
# ls NOVOPlasty/*.pl

# # Seed: the Vorticella cytochrome b contig
# module load seqkit
# seqkit grep -p k141_466 ../megahit_out2/final.contigs.fa | seqkit seq -i > seed.fa
# seqkit stats seed.fa

