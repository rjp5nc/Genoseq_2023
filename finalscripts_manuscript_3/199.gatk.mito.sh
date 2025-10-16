


module load bcftools

cd /scratch/rjp5nc/UK2022_2024/mito_vcf/merged_gvcf

bcftools query -l /scratch/rjp5nc/UK2022_2024/mito_vcf/merged_gvcf/usdobtusa_mito_combined.g.vcf.gz > samples_mito_USobtusa_snapp.txt
