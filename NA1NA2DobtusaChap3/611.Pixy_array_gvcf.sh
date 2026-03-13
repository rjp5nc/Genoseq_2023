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






gawk -v OFS="\t" '
BEGINFILE {
  if (ARGIND == 1) FS="\t";
  else if (ARGIND == 2 || ARGIND == 3) FS=",";
  else if (ARGIND == 4) FS="[ \t]+";   # <-- robust: tabs OR spaces
}

# 1) pops: sample \t pop
ARGIND==1 { pop[$1] = $2; next }

# 2) genomic_types_v2.csv: pick 1 rep per Group
ARGIND==2 {
  gsub(/"/, "", $0)
  if (FNR == 1) next
  clone = $3
  grp   = $4
  if (!(seen_grp[grp]++)) {
    rep[grp] = clone
    order[++n] = grp
  }
  next
}

# 3) mito_types.csv
ARGIND==3 {
  gsub(/"/, "", $0)
  if (FNR == 1) next
  mito[$1] = $2
  next
}

# 4) OO groups: always include these samples
ARGIND==4 {
  oo_sample = $1
  oo_grp    = $2
  forced_s[++nf] = oo_sample
  forced_g[oo_sample] = oo_grp
  next
}

END {
  for (i=1; i<=n; i++) {
    grp = order[i]
    s   = rep[grp]
    printed[s] = 1

    p = (s in pop  ? pop[s]  : "UnknownPop")
    m = (s in mito ? mito[s] : "UnknownMito")
    print s, p "_" grp "_" m
  }

  for (j=1; j<=nf; j++) {
    s = forced_s[j]
    if (s in printed) continue
    grp = forced_g[s]

    p = (s in pop  ? pop[s]  : "UnknownPop")
    m = (s in mito ? mito[s] : "UnknownMito")
    print s, p "_" grp "_" m
  }
}
' \
/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/pops_pixy_ready.txt \
/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/genomic_types_v2.csv \
/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/mito_types.csv \
/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/oo_groups.tsv \
> /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/one_per_superclone_plus_OO_with_mitotype.txt

grep -v "Unknown" \
  /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/one_per_superclone_plus_OO_with_mitotype.txt \
> /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/one_per_superclone_plus_OO_with_mitotype.noUnknown.txt

# Directory with chromosome VCFs
VCF_DIR=/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/usdobtusa_gvcf_10bp/usobtusa_vcf_snps
OUT_DIR=/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/results_pixy_chr_2
mkdir -p $OUT_DIR

# Get list of chromosome VCFs
VCFS=($VCF_DIR/*_filtsnps10bpindels_snps.vcf.gz)

# Select VCF for this array task
CHR_VCF=${VCFS[$SLURM_ARRAY_TASK_ID]}
CHR_NAME=$(basename $CHR_VCF _filtsnps10bpindels_snps.vcf.gz)

echo "Running Pixy on $CHR_NAME"

pixy --stats pi fst dxy \
    --vcf $CHR_VCF \
    --populations /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/one_per_superclone_plus_OO_with_mitotype.noUnknown.txt \
    --window_size 10000 \
    --n_cores 20 \
    --output_folder $OUT_DIR \
    --output_prefix pixy_${CHR_NAME}



# cat /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/results_pixy_chr/*_fst.txt > /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/usobtusa_fst.txt
# cat /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/results_pixy_chr/*_dxy.txt > /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/usobtusa_dxy.txt
# cat /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/results_pixy_chr/*_pi.txt > /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/usobtusa_pi.txt