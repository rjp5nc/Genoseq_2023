#!/usr/bin/env bash

#SBATCH -J BEAST # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 6-0:00:00 ### 15 seconds
#SBATCH --mem 60G
#SBATCH -o /scratch/rjp5nc/erroroutputs/beast.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/erroroutputs/beast.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab
#SBATCH --mail-type=END               # Send email at job completion
#SBATCH --mail-user=rjp5nc@virginia.edu    # Email address for notifications

# This script will run snpEFF to annotate a VCF.

# Load Modules
module load tabix/0.2.6
module load gcc/9.2.0 htslib/1.10.2

# Working folder is core folder where this pipeline is being run.
wd="/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/Daphnia_obtusa_genes/"

# Combined VCF name
vcf="/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/trimmed10bp_filtered_Gilmer_12_renamed.vcf.gz"

# Output annotated VCF name
out_vcf="output.ann.vcf"

# java parameters
JAVAMEM=50G

# Move to working directory
cd ${wd}

# Run snpEFF on raw vcf
java -Xmx${JAVAMEM} -jar \
/home/csm6hg/SNPEFF/snpEff.jar ann \
dpgenome \
$vcf \
-o vcf \
-t \
-s ${wd}/SNPeff.summary > \
$out_vcf

# BGzip and tab index annotated vcf
bgzip ${wd}/${out_vcf}
tabix -p vcf ${wd}/${out_vcf}.gz

# Finish
echo "Complete" $(date)







~/miniconda3/envs/snpeff_env/bin/snpEff \
-c /scratch/rjp5nc/snpEff_notworking/snpEff.config \
US_obtusa \
"$vcf" \
-o vcf \
-s "${wd}/SNPeff.summary" > "$out_vcf"