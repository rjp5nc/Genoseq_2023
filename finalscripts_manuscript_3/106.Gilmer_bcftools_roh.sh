#!/usr/bin/env bash
#
#SBATCH -J bcftools_roh # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-48:00:00  ### 48 hours
#SBATCH --mem 50G
#SBATCH -o /scratch/rjp5nc/err/boh.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/err/boh.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab



# conda create -n rohviz -c conda-forge -c bioconda bcftools perl -y
# conda activate rohviz

# curl -L -o "$CONDA_PREFIX/bin/roh-viz" \
#   https://raw.githubusercontent.com/samtools/bcftools/refs/heads/develop/misc/roh-viz
# chmod +x "$CONDA_PREFIX/bin/roh-viz"

# # sanity check
# which roh-viz



conda init bash
source ~/.bashrc
conda activate rohviz


# module purge
# module spider bcftools         # see available versions
# module load bcftools/1.17

module purge

module load bcftools/1.17

dir=/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv
vcf=trimmed10bp_filtered_Gilmer.vcf.gz
txt=bcftools_roh_gilmer_out.txt

cd $dir

bcftools roh -G30 trimmed10bp_filtered_Gilmer.vcf.gz -o bcftools_roh_gilmer_out.txt
$CONDA_PREFIX/bin/roh-viz -i ${dir}/$txt -v ${dir}/$vcf -o ${dir}/rmme.html
