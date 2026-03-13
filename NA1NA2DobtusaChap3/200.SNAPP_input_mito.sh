#!/usr/bin/env bash
#
#SBATCH -J Snapp_input # A single job name for the array
#SBATCH --ntasks-per-node=15 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 6-0:00  ### 48 hours
#SBATCH --mem 200G
#SBATCH -o /scratch/rjp5nc/err/trim.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/err/trim.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab
#!/bin/bash
# Script used to prepare mitochondrial VCF and run SNAPP for phylogenetic splits


# Load modules
module load bcftools
module load ruby




# Working directory
wd="/scratch/rjp5nc/snapp5"
cd $wd


bcftools query -f '[%SAMPLE\t%DP\n]' /scratch/rjp5nc/UK2022_2024/mito_vcf/merged_gvcf/usdobtusa_mito_combined.g.vcf.gz \
> usdobtusa_mito_DP_per_sample.txt

awk '
{
  sum[$1] += $2;       # add DP to sample sum
  count[$1] += 1;      # count sites for sample
}
END {
  for (s in sum) {
    print s, sum[s]/count[s]   # average DP per sample
  }
}' usdobtusa_mito_DP_per_sample.txt > usdobtusa_mito_avg_DP_per_sample.txt

# awk '$2>10 {print $1}' usdobtusa_mito_avg_DP_per_sample.txt > high_DP_samples.txt


# bcftools view -S high_DP_samples.txt \
#   /scratch/rjp5nc/UK2022_2024/mito_vcf/merged_gvcf/usdobtusa_mito_combined.g.vcf.gz \
#   -Oz -o usdobtusa_mito_highDP_samples.g.vcf.gz
# bcftools index usdobtusa_mito_highDP_samples.g.vcf.gz
# tabix -p vcf usdobtusa_mito_highDP_samples.g.vcf.gz


# filter VCF by samples and bgzip output
bcftools view -S /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/final_vcf_filter_one_of_each.txt -Oz -o usdobtusa_mito_highDP_samples2.g.vcf.gz /scratch/rjp5nc/UK2022_2024/mito_vcf/merged_gvcf/usdobtusa_mito_combined.g.vcf.gz

# index the new VCF
bcftools index -f usdobtusa_mito_highDP_samples2.g.vcf.gz
tabix -p vcf usdobtusa_mito_highDP_samples2.g.vcf.gz


# Mitochondrial VCF
vcf=${wd}/usdobtusa_mito_highDP_samples2.g.vcf.gz


gatk GenotypeGVCFs \
  -R /scratch/rjp5nc/Reference_genomes/mito_reference/usdobtusa_mito.fasta \
  -V ${vcf} \
  -O usdobtusa_mito_joint.vcf.gz

bcftools view \
  --threads 15 \
  -m2 -M2 \
  -v snps \
  -Oz \
  -o usdobtusa_mito_biallelic.vcf.gz \
  usdobtusa_mito_joint.vcf.gz

bcftools index usdobtusa_mito_biallelic.vcf.gz

# Beast2 directory
beast2="/scratch/rjp5nc/beast/beast/bin/beast"

# Sample/group information
# (You can adjust these if you want to use the same grouping as nuclear data)
IN=/scratch/rjp5nc/UK2022_2024/mito_vcf/merged_gvcf/samples_mito_USobtusa_snapp.txt
OUTTMP=/scratch/rjp5nc/UK2022_2024/mito_vcf/usdobtusa_mito_clone.tmp.txt

# Clean VCF (remove ALT="*")
bcftools view \
  -e 'ALT="*"' \
  -m2 -M2 \
  -Ov \
  -o ${wd}/usdobtusa.mito.biallelic.clean.vcf \
  ${wd}/usdobtusa_mito_biallelic.vcf.gz





bcftools query -l /scratch/rjp5nc/snapp5/usdobtusa.mito.biallelic.clean.vcf \
  > /scratch/rjp5nc/snapp5/samples_mito_USobtusa_snapp_clean.txt

echo -n "monophyletic NA " > snapp_constraints.txt
sed 's/$/_clone/' samples_mito_USobtusa_snapp_clean.txt | paste -sd, - >> snapp_constraints.txt

awk '{print $1 "_clone", $1}' /scratch/rjp5nc/snapp5/samples_mito_USobtusa_snapp_clean.txt > samples_mito_USobtusa_snapp_clean_2col.txt



usdobtusa_mito_avg_DP_per_sample.txt






ruby ${wd}/snapp_prep.rb \
  -v ${wd}/usdobtusa.mito.biallelic.clean.vcf \
  -t ${wd}/samples_mito_USobtusa_snapp_clean_2col.txt \
  -c snapp_constraints.txt \
  -x ${wd}/snapp.mito.xml \
  -o ${wd}/snapp.mito \
  -m 1000 \
  -l 10000


#INFO: Retained 156 bi-allelic sites. at dp > 10

# Optional: Run Beast2 (uncomment if ready)
# ${beast2} ${wd}/snapp.mito.xml