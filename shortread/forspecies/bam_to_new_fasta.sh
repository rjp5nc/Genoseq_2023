#!/bin/bash

module load samtools

mkdir -p "$output_folder"

for bam in /scratch/rjp5nc/UK2022_2024/allshortreads/mitosortedbams/*.bam; do
    base_name=$(basename "$bam" .bam)
    samtools fastq "$bam" -o "/scratch/rjp5nc/UK2022_2024/allshortreads/mitofastas/${base_name}_allreads.fastq"
done


for fastq in /scratch/rjp5nc/UK2022_2024/allshortreads/mitofastas/*.fastq; do
    base_name=$(basename "$fastq" .fastq)
    bwa mem -t 10 -R "@RG\tID:${base_name}\tSM:${base_name}\tPL:ILLUMINA" \
        /scratch/rjp5nc/UK2022_2024/allshortreads/mtdna_D8_119.fa "$fastq" | \
        samtools view -bS -o "/scratch/rjp5nc/UK2022_2024/allshortreads/alignedbamsnew/${base_name}.bam"
done


for bam in /scratch/rjp5nc/UK2022_2024/allshortreads/alignedbamsnew/*.sorted_allreads.bam; do
    base_name=$(basename "$bam" .sorted_allreads.bam)
  
  # Define the output filtered BAM file
  output_file="/scratch/rjp5nc/UK2022_2024/allshortreads/filteredbams/${base_name}.filtered.bam"
  
  # Run samtools to filter out unmapped reads (-F 4)
  samtools view -b -F 4 "$bam" > "$output_file"
  
  echo "Filtered BAM saved to: $output_file"
done




BAM_DIR="/scratch/rjp5nc/UK2022_2024/allshortreads/filteredbams"

# Loop through all BAM files in the directory
for bam_file in "$BAM_DIR"/*.bam; do
    echo "Indexing $bam_file"
    samtools index "$bam_file"
done





bam_list=$(ls /scratch/rjp5nc/UK2022_2024/allshortreads/filteredbams/*.filtered.bam | sort)

# Generate VCF with VarScan
samtools mpileup \
    --fasta-ref /scratch/rjp5nc/UK2022_2024/allshortreads/mtdna_D8_119.fa \
    $bam_list | \

java -jar $EBROOTVARSCAN/VarScan.v2.4.4.jar mpileup2snp \
    /dev/stdin \
    --min-coverage 4 \
    --min-var-freq 0.001 \
    --output-vcf > /scratch/rjp5nc/UK2022_2024/allshortreads/mitoVCF2.vcf

# Extract sample names in correct order
echo "$bam_list" | sed 's#.*/##; s/.bam//' > sample_names.txt

# Reheader the VCF with the correct sample order
bcftools reheader \
    --samples sample_names.txt \
    /scratch/rjp5nc/UK2022_2024/allshortreads/mitoVCF2.vcf > \
    /scratch/rjp5nc/UK2022_2024/allshortreads/mitoVCF_renamed2.vcf