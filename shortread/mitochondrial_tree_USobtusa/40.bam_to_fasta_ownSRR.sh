#!/usr/bin/env bash

#SBATCH -J bamtofasta # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 2:00:00 ### 15 seconds
#SBATCH --mem 50G
#SBATCH -o /scratch/rjp5nc/erroroutputs/nFlo_1.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/erroroutputs/nFlo_1.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab
#SBATCH --array=1-22
#SBATCH --mail-type=END               # Send email at job completion
#SBATCH --mail-user=rjp5nc@virginia.edu    # Email address for notifications


module load bcftools
module load htslib 
module load gsl
module load samtools

#SLURM_ARRAY_TASK_ID=1

cd /scratch/rjp5nc/UK2022_2024/final_mitobam_rg2
outfq=/scratch/rjp5nc/UK2022_2024/final_mitobam_rg2
outfq2=/scratch/rjp5nc/UK2022_2024/consensusmito

CSV_FILE="/scratch/rjp5nc/UK2022_2024/forref2_mito.csv"

line=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" ${CSV_FILE})

# Extract fields (assuming CSV format: sample_id,reference_path)
samp=$(echo "$line" | cut -d',' -f4)
ref_path=$(echo "$line" | cut -d',' -f5)
ref_path=$(echo "${ref_path}" | tr -d '\r')
# Call snps
bcftools mpileup -q 30 -Q 20 -Ou -f ${ref_path} ${outfq}/${samp}finalmap_RG.bam | \
bcftools call -mv -V indels | \
bcftools filter -i 'QUAL>20 && INFO/DP>=20' -Oz -o ${outfq}/${samp}.filt.mito.vcf.gz

  # Index VCF
  tabix -p vcf ${outfq}/${samp}.filt.mito.vcf.gz

  # apply variants to create consensus sequence
  cat ${ref_path} | \
  bcftools consensus ${outfq}/${samp}.filt.mito.vcf.gz \
  --sample ${samp} > \
  ${outfq2}/${samp}.filt.consensus.mito.fa

  # Change name of header
  name=$( echo ${samp}.filt.consensus.mito.fa | tr "." '\t' | cut -f1 | tr "/" '\t' | cut -f1 )
  echo ${name}

  # Extracts chromosome via bed file
  sed -i "s/^>*"mtdna"/>"mtdna."${name}/g" \
  ${outfq2}/${samp}.filt.consensus.mito.fa


  # Index fasta
  samtools faidx ${outfq2}/${samp}.filt.consensus.mito.fa