#!/usr/bin/env bash
#
#SBATCH -J vort_mito_asm
#SBATCH --ntasks-per-node=8
#SBATCH -N 1
#SBATCH -t 24:00:00
#SBATCH --mem 120G
#SBATCH -o /scratch/rjp5nc/err/vort_mito_asm.%A_%a.out
#SBATCH -e /scratch/rjp5nc/err/vort_mito_asm.%A_%a.err
#SBATCH -p standard
#SBATCH --account berglandlab_standard
#SBATCH --mail-type=END
#SBATCH --mail-user=rjp5nc@virginia.edu








# cd /scratch/rjp5nc/vorticella
# mkdir -p vort_mito_asm && cd vort_mito_asm

# # Get NOVOPlasty
# git clone https://github.com/ndierckx/NOVOPlasty.git
# ls NOVOPlasty/*.pl

# # Seed: the Vorticella cytochrome b contig
# module load seqkit
# seqkit grep -p k141_466 ../megahit_out2/final.contigs.fa | seqkit seq -i > seed.fa
# seqkit stats seed.fa



# cd /scratch/rjp5nc/vorticella/vort_mito_asm
# cat > config.txt << 'EOF'
# Project:
# -----------------------
# Project name          = Vort_mito
# Type                  = mito
# Genome Range          = 25000-70000
# K-mer                 = 33
# Max memory            = 100
# Extended log          = 0
# Save assembled reads  = no
# Seed Input            = /scratch/rjp5nc/vorticella/vort_mito_asm/seed.fa
# Extend seed directly  = no
# Reference sequence    = 
# Variance detection    = 
# Chloroplast sequence  = 

# Dataset 1:
# -----------------------
# Read Length           = 150
# Insert size           = 300
# Platform              = illumina
# Single/Paired         = PE
# Combined reads        = 
# Forward reads         = /scratch/rjp5nc/vorticella/nonDaphnia2_R1.fq.gz
# Reverse reads         = /scratch/rjp5nc/vorticella/nonDaphnia2_R2.fq.gz
# Store Hash            = 

# Heteroplasmy:
# -----------------------
# MAF                   = 
# HP exclude list       = 
# PCR-free              = 

# Optional:
# -----------------------
# Insert size auto      = yes
# Use Quality Scores    = no
# Reduce ambigious N's  = 
# Output path           = /scratch/rjp5nc/vorticella/vort_mito_asm/out/
# EOF

# mkdir -p out
# head -12 config.txt
# ls -lh /scratch/rjp5nc/vorticella/nonDaphnia2_R1.fq.gz /scratch/rjp5nc/vorticella/nonDaphnia2_R2.fq.gz




set -euo pipefail
cd /scratch/rjp5nc/vorticella/vort_mito_asm
mkdir -p out

novo=$(ls NOVOPlasty/NOVOPlasty*.pl | head -1)
echo "Using $novo"
perl "$novo" -c config.txt

echo "Done: $(date)"
ls -lh out/
