#!/usr/bin/env bash
#
#SBATCH -J pixy # A single job name for the array
#SBATCH --ntasks-per-node=10 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-72:00  ### 48 hours
#SBATCH --mem 10G
#SBATCH -o /scratch/rjp5nc/err/pixy.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/err/pixy.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab


#ijob -A berglandlab -c10 -p standard --mem=40G

#conda create -n pixy -c conda-forge -c bioconda pixy -y
# activate the environment
conda activate pixy

#conda install --yes -c conda-forge pixy
#conda install -c conda-forge pixy=2.0.0.beta11

#conda install -c bioconda htslib
#conda install -c bioconda vcftools


module load bcftools
cd /scratch/rjp5nc/UK2022_2024/mito_vcf/

awk -F',' 'NR>1 {print $1 "\t" $6}' /project/berglandlab/Robert/UKSequencing2022_2024/old_stuff/2022_2024seqmetadata20250131.csv > pops.txt

awk 'NF==2' pops.txt > pops_filtered.txt


bcftools query -l  /scratch/rjp5nc/UK2022_2024/mito_vcf/merged_vcfs/usdobtusa_mito_genotyped.vcf.gz > vcf_samples.txt

grep -Ff vcf_samples.txt pops_filtered.txt > pops_both.txt

grep -Ff <(bcftools query -l /scratch/rjp5nc/UK2022_2024/mito_vcf/merged_vcfs/usdobtusa_mito_genotyped.vcf.gz) pops_both.txt > pops_fixed.txt


VCF="/scratch/rjp5nc/UK2022_2024/mito_vcf/merged_vcfs/usdobtusa_mito_genotyped.vcf.gz"
POPS="/scratch/rjp5nc/UK2022_2024/mito_vcf/pops_fixed.txt"

# --- Output directory ---
OUTDIR="/scratch/rjp5nc/UK2022_2024/mito_vcf/results_pixy"
mkdir -p ${OUTDIR}

# --- Parameters ---
WINDOW=500   # window size in bp, adjust as needed

cd /scratch/rjp5nc/UK2022_2024/mito_vcf
tabix -p vcf /scratch/rjp5nc/UK2022_2024/mito_vcf/merged_vcfs/usdobtusa_mito_genotyped.vcf.gz
chmod u+rwx /scratch/rjp5nc/UK2022_2024/mito_vcf/results_pixy

cut -f1 pops_fixed.txt > pops_samples.txt


bcftools query -l /scratch/rjp5nc/UK2022_2024/mito_vcf/merged_vcfs/usdobtusa_mito_genotyped.vcf.gz > vcf_samples.txt

# Keep only matching samples
grep -F -x -f vcf_samples.txt pops_samples.txt > pops_samples_filtered.txt

bcftools view -S pops_samples_filtered.txt -Oz -o usdobtusa_mito_genotyped_subset.vcf.gz /scratch/rjp5nc/UK2022_2024/allsites_mito/usdobtusa_mito_allsites.vcf.gz

grep -F -w -f pops_samples_filtered.txt pops_fixed.txt > pops_fixed_filtered.txt

tabix -p vcf usdobtusa_mito_genotyped_subset.vcf.gz

chmod u+rwx /scratch/rjp5nc/UK2022_2024/mito_vcf/results_pixy

bcftools query -f '[%DP\t]\n' usdobtusa_mito_genotyped_subset.vcf.gz | \
awk '
{for(i=1;i<=NF;i++){sum[i]+=$i; count[i]++}}
END{
    for(i=1;i<=NF;i++){
        mean=sum[i]/count[i];
        print i, mean
    }
}' > sample_mean_depth.txt

bcftools query -l usdobtusa_mito_genotyped_subset.vcf.gz > sample_names.txt
paste sample_names.txt sample_mean_depth.txt > sample_depths_mapped.txt

awk '$3<2 {print $1}' sample_depths_mapped.txt > low_coverage_samples.txt

bcftools view -s ^$(paste -sd, low_coverage_samples.txt) \
    -Oz -o usdobtusa_mito_filtered.vcf.gz \
    usdobtusa_mito_genotyped_subset.vcf.gz

bcftools index -t usdobtusa_mito_filtered.vcf.gz
bcftools query -l usdobtusa_mito_filtered.vcf.gz > filtered_samples.txt
grep -Ff filtered_samples.txt pops_fixed_filtered.txt > pops_filtered_for_pixy.txt

bcftools query -l usdobtusa_mito_filtered.vcf.gz > vcf_samples.txt
awk -F'\t' 'NR==FNR {vcf[$1]; next} $1 in vcf' vcf_samples.txt pops_fixed_filtered.txt > pops_pixy_ready.txt


pixy --stats pi fst dxy \
--vcf usdobtusa_mito_filtered.vcf.gz \
--populations pops_pixy_ready.txt \
--window_size 100 \
--n_cores 4 \
--output_folder /scratch/rjp5nc/UK2022_2024/mito_vcf/results_pixy/ \
--bypass_invariant_check \
--output_prefix pixy

# --- Done ---
echo "✅ Pixy run finished. Outputs in ${OUTDIR}/"


#awk '{print $1 > $2".txt"}' pops_fixed_filtered.txt

#for popfile in P58.txt P62.txt P63.txt P66.txt Gilmer.txt; do
#    popname=$(basename "$popfile" .txt)
#    vcftools --gzvcf usdobtusa_mito_genotyped_subset.vcf.gz \
#             --keep "$popfile" \
#             --site-pi \
#             --out "${popname}_pi"
#done
