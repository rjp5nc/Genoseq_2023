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

vcf=/scratch/rjp5nc/UK2022_2024/allsites_mito/usdobtusa_mito_allsites.haploid.vcf.gz
bcftools query -l "$vcf" > /scratch/rjp5nc/UK2022_2024/allsites_mito/mito_samples.txt

contig=CM028013.1
len=$(bcftools view -h "$vcf" | awk -v c=$contig '$0 ~ "##contig=<ID="c"," {match($0,/length=([0-9]+)/,m); if(m[1]) print m[1];}')
printf "%s\t0\t%s\n" "$contig" "$len" > /scratch/rjp5nc/UK2022_2024/allsites_mito/mito.whole.bed



awk -F',' 'NR>1 {print $1 "\t" $1}' /project/berglandlab/Robert/UKSequencing2022_2024/old_stuff/2022_2024seqmetadata20250131.csv > pops.txt

#awk -F',' 'NR>1 {print $1 "\t" $6}' /project/berglandlab/Robert/UKSequencing2022_2024/old_stuff/2022_2024seqmetadata20250131.csv > pops.txt

awk 'NF==2' pops.txt > pops_filtered.txt

# /scratch/rjp5nc/UK2022_2024/daphnia_phylo/Sample_ID_Species_merged_20250627.csv
echo -e "SRR5012393\tSRR5012393\nSRR5012394\tSRR5012394\nSRR5012770\tSRR5012770\nSRR5012396\tSRR5012396\nSRR5012773\tSRR5012773\nSRR5012771\tSRR5012771" \
>> /scratch/rjp5nc/UK2022_2024/mito_vcf/pops_filtered.txt


$CONDA_PREFIX/bin/bcftools query -l  $vcf > vcf_samples.txt

grep -Ff vcf_samples.txt pops_filtered.txt > pops_both.txt

awk '$2 != "PBO66"' pops_both.txt > pops_both2.txt

grep -wFf <($CONDA_PREFIX/bin/bcftools query -l $vcf) \
  pops_both2.txt > pops_fixed.txt




POPS=/scratch/rjp5nc/UK2022_2024/mito_vcf/pops_fixed.txt

OUT=/scratch/rjp5nc/UK2022_2024/allsites_mito/usdobtusa_mito_subset.vcf.gz
awk '{print $1}' "$POPS" | sort | uniq > /scratch/rjp5nc/UK2022_2024/allsites_mito/pixy_keep_samples.txt



bcftools view \
  -S /scratch/rjp5nc/UK2022_2024/allsites_mito/pixy_keep_samples.txt \
  -Oz -o "$OUT" "$vcf"
  
bcftools index -t "$OUT"

conda init bash
source ~/.bashrc  # or restart your shell
conda activate pixy
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH

pixy --stats pi dxy \
     --vcf "$OUT" \
     --populations /scratch/rjp5nc/UK2022_2024/mito_vcf/pops_fixed.txt \
     --bed_file /scratch/rjp5nc/UK2022_2024/allsites_mito/mito.whole.bed \
     --n_cores 8 \
         --output_folder /scratch/rjp5nc/UK2022_2024/allsites_mito/pixy_mito_whole \
    --output_prefix pixy_14601_incSRR


