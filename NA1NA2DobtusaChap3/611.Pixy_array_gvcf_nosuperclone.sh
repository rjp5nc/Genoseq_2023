#!/usr/bin/env bash
#
#SBATCH -J pixy # A single job name for the array
#SBATCH --ntasks-per-node=10 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-92:00  ### 48 hours
#SBATCH --mem 10G
#SBATCH -o /scratch/rjp5nc/err/pixy.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/err/pixy.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab
#SBATCH --array=1-135%70

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

# gawk -v OFS="\t" '
# BEGINFILE {
#   if (ARGIND == 1) FS="\t";
#   else if (ARGIND == 2 || ARGIND == 3) FS=",";
#   else if (ARGIND == 4) FS="[ \t]+";
# }

# # 1) pops: sample \t pool
# ARGIND==1 { pool[$1] = $2; next }

# # 2) genomic_types_v2.csv: pick 1 rep per superclone Group
# ARGIND==2 {
#   gsub(/"/, "", $0)
#   if (FNR == 1) next
#   clone = $3
#   grp   = $4
#   if (!(seen_grp[grp]++)) {
#     rep[grp] = clone
#     order[++n] = grp
#   }
#   next
# }

# # 3) mito_types.csv: sampleA,Group
# ARGIND==3 {
#   gsub(/"/, "", $0)
#   if (FNR == 1) next
#   mito[$1] = $2
#   next
# }

# # 4) OO groups: always include these samples (no need to store grp in output)
# ARGIND==4 {
#   forced[++nf] = $1
#   next
# }

# END {
#   # reps from genomic groups
#   for (i=1; i<=n; i++) {
#     s = rep[order[i]]
#     printed[s] = 1

#     p = (s in pool ? pool[s] : "UnknownPool")
#     m = (s in mito ? mito[s] : "UnknownMito")

#     if (p ~ /Unknown/ || m ~ /Unknown/) continue
#     print s, p "_" m
#   }

#   # forced OO samples (skip duplicates)
#   for (j=1; j<=nf; j++) {
#     s = forced[j]
#     if (s in printed) continue

#     p = (s in pool ? pool[s] : "UnknownPool")
#     m = (s in mito ? mito[s] : "UnknownMito")

#     if (p ~ /Unknown/ || m ~ /Unknown/) continue
#     print s, p "_" m
#   }
# }
# ' \
# pops_pixy_ready.txt \
# genomic_types_v2.csv \
# /scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/mito_types.csv \
# oo_groups.tsv \
# > one_per_superclone_pool_mito.txt


# Directory with chromosome VCFs
VCF_DIR=/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/usdobtusa_gvcf_10bp/usobtusa_vcf_snps
OUT_DIR=/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/results_pixy_chr_3
mkdir -p $OUT_DIR

# Get list of chromosome VCFs
VCFS=($VCF_DIR/*_filtsnps10bpindels_snps.vcf.gz)

# Select VCF for this array task
CHR_VCF=${VCFS[$SLURM_ARRAY_TASK_ID]}
CHR_NAME=$(basename $CHR_VCF _filtsnps10bpindels_snps.vcf.gz)

echo "Running Pixy on $CHR_NAME"

pixy --stats pi fst dxy \
    --vcf $CHR_VCF \
    --populations /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/one_per_superclone_pool_mito.txt \
    --window_size 10000 \
    --n_cores 20 \
    --output_folder $OUT_DIR \
    --output_prefix pixy_${CHR_NAME}



# cat /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/results_pixy_chr_2/*_fst.txt > /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/usobtusa_fst_mitosup.txt
# cat /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/results_pixy_chr_2/*_dxy.txt > /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/usobtusa_dxy_mitosup.txt
# cat /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/results_pixy_chr_2/*_pi.txt > /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/usobtusa_pi_mitosup.txt