#!/usr/bin/env bash
#
#SBATCH -J pixy # A single job name for the array
#SBATCH --ntasks-per-node=20 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-92:00  ### 48 hours
#SBATCH --mem 80G
#SBATCH -o /scratch/rjp5nc/err/pixy.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/err/pixy.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

#ijob -A berglandlab -c10 -p standard --mem=40G

#conda create -n pixy -c conda-forge -c bioconda pixy -y
# activate the environment
conda init bash
source ~/.bashrc  # or restart your shell
conda activate pixy
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH

#conda install --yes -c conda-forge pixy
#conda install -c conda-forge pixy=2.0.0.beta11

#conda install -c bioconda htslib
#conda install -c bioconda vcftools
#conda install -c bioconda bcftools

cd /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv

file /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/trimmed10bp_allsites_Repeatmasked_usobtusa.filtered_bgz2.vcf.gz


awk -F',' 'NR>1 {print $1 "\t" $6}' /project/berglandlab/Robert/UKSequencing2022_2024/old_stuff/2022_2024seqmetadata20250131.csv > pops.txt

awk 'NF==2' pops.txt > pops_filtered.txt

$CONDA_PREFIX/bin/bcftools query -l  /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/trimmed10bp_allsites_Repeatmasked_usobtusa.filtered_bgz2.vcf.gz > vcf_samples.txt

grep -Ff vcf_samples.txt pops_filtered.txt > pops_both.txt

awk '$2 != "PBO66"' pops_both.txt > pops_both2.txt

grep -Ff <($CONDA_PREFIX/bin/bcftools query -l trimmed10bp_allsites_Repeatmasked_usobtusa.filtered_bgz2.vcf.gz) pops_both2.txt > pops_fixed.txt


VCF="/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/trimmed10bp_allsites_Repeatmasked_usobtusa.filtered_bgz2.vcf.gz"
POPS="/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/pops_fixed.txt"

# --- Output directory ---
OUTDIR="/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/results_pixy10000"
mkdir -p ${OUTDIR}

# --- Parameters ---
WINDOW=10000   # window size in bp, adjust as needed

chmod u+rwx /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/results_pixy10000

