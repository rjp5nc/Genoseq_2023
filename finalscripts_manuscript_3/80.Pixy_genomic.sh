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

file trimmed10bp_allsites_usobtusa.vcf.gz


awk -F',' 'NR>1 {print $1 "\t" $6}' /project/berglandlab/Robert/UKSequencing2022_2024/old_stuff/2022_2024seqmetadata20250131.csv > pops.txt

awk 'NF==2' pops.txt > pops_filtered.txt

$CONDA_PREFIX/bin/bcftools query -l  trimmed10bp_allsites_usobtusa.vcf.gz > vcf_samples.txt

grep -Ff vcf_samples.txt pops_filtered.txt > pops_both.txt

awk '$2 != "PBO66"' pops_both.txt > pops_both2.txt

grep -Ff <($CONDA_PREFIX/bin/bcftools query -l trimmed10bp_allsites_usobtusa.vcf.gz) pops_both2.txt > pops_fixed.txt


VCF="/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/trimmed10bp_allsites_usobtusa.vcf.gz"
POPS="/scratch/rjp5nc/UK2022_2024/mito_vcf/pops_fixed.txt"

# --- Output directory ---
OUTDIR="/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/results_pixy"
mkdir -p ${OUTDIR}

# --- Parameters ---
WINDOW=10000   # window size in bp, adjust as needed

chmod u+rwx /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/results_pixy

cut -f1 pops_fixed.txt > pops_samples.txt

# Keep only matching samples
grep -F -x -f vcf_samples.txt pops_samples.txt > pops_samples_filtered.txt

#bcftools view -S pops_samples_filtered.txt -Oz \
#    --threads 16 \
#    -o trimmed10bp_allsites_usobtusa.filtered_bgz.vcf.gz \
#    trimmed10bp_allsites_usobtusa.bgz.vcf.gz

tabix -p vcf trimmed10bp_allsites_usobtusa.filtered_bgz.vcf.gz

$CONDA_PREFIX/bin/bcftools index -t trimmed10bp_allsites_usobtusa.filtered_bgz.vcf.gz
$CONDA_PREFIX/bin/bcftools query -l trimmed10bp_allsites_usobtusa.filtered_bgz.vcf.gz > filtered_samples.txt
grep -Ff filtered_samples.txt pops_fixed.txt > pops_filtered_for_pixy.txt

$CONDA_PREFIX/bin/bcftools query -l trimmed10bp_allsites_usobtusa.filtered_bgz.vcf.gz > vcf_samples2.txt
awk -F'\t' 'NR==FNR {vcf[$1]; next} $1 in vcf' vcf_samples2.txt pops_filtered_for_pixy.txt > pops_pixy_ready.txt


