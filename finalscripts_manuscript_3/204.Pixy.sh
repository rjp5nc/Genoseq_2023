#!/usr/bin/env bash
#
#SBATCH -J pixy # A single job name for the array
#SBATCH --ntasks-per-node=20 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-92:00  ### 48 hours
#SBATCH --mem 50G
#SBATCH -o /scratch/rjp5nc/err/pixy.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/err/pixy.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

#ls $VCF_DIR/*_filtsnps10bpindels_snps.vcf.gz | wc -l
###SBATCH --array=1-347%50

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

cd /scratch/rjp5nc/UK2022_2024/mito_vcf

# Directory with chromosome VCFs

VCF="/scratch/rjp5nc/UK2022_2024/allsites_mito/usdobtusa_mito_allsites_firstcol.vcf.gz"
POPS="/scratch/rjp5nc/UK2022_2024/allsites_mito/pops_fixed.txt"

# --- Output directory ---
OUT_DIR="/scratch/rjp5nc/UK2022_2024/allsites_mito/results_pixymito"

pixy --stats pi fst dxy \
    --vcf $VCF \
    --populations /scratch/rjp5nc/UK2022_2024/allsites_mito/pops_fixed.txt \
    --window_size 14601 \
    --n_cores 20 \
    --output_folder $OUT_DIR \
    --bypass_invariant_check \
    --output_prefix pixy_14601_


# cat /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/results_pixy_chr/*_fst.txt > /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/usobtusa_fst.txt
# cat /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/results_pixy_chr/*_dxy.txt > /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/usobtusa_dxy.txt
# cat /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/results_pixy_chr/*_pi.txt > /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/usobtusa_pi.txt