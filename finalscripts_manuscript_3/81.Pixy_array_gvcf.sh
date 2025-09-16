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
#SBATCH --array=1-347

#ls $VCF_DIR/*_filtsnps10bpindels_snps.vcf.gz | wc -l

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


# Directory with chromosome VCFs
VCF_DIR=/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/usdobtusa_gvcf_10bp/usobtusa_vcf_snps/
OUT_DIR=/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/results_pixy_chr
mkdir -p $OUT_DIR

# Get list of chromosome VCFs
VCFS=($VCF_DIR/*_filtsnps10bpindels_snps.vcf.gz)

# Select VCF for this array task
CHR_VCF=${VCFS[$SLURM_ARRAY_TASK_ID]}
CHR_NAME=$(basename $CHR_VCF _filtsnps10bpindels_snps.vcf.gz)

echo "Running Pixy on $CHR_NAME"

pixy --stats pi fst dxy \
    --vcf $CHR_VCF \
    --populations /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/pops_pixy_ready.txt \
    --window_size 100 \
    --n_cores 20 \
    --output_folder $OUT_DIR \
    --output_prefix pixy_${CHR_NAME}