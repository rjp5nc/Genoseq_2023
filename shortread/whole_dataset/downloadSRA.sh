#!/usr/bin/env bash
#
#SBATCH -J vcf2gds # A single job name for the array
#SBATCH --ntasks-per-node=5 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-48:00  ### 48 hours
#SBATCH --mem 20G
#SBATCH -o /scratch/rjp5nc/downloadsra.out # Standard output
#SBATCH -e /scratch/rjp5nc/downloadsra.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab



#du -sh /scratch/rjp5nc
#ls /project/berglandlab/daphnia_genus/short_read/murray_data/all_gvcf/Euro_gvcfs > /scratch/rjp5nc/Euro_gvcfs.txt
#ls /project/berglandlab/daphnia_genus/short_read/murray_data/all_gvcf/gvcf/Scaffold_9201_HRSCAF_10758 > /scratch/rjp5nc/Scaffold_9201_HRSCAF_10758.txt
#ls /project/berglandlab/daphnia_genus/short_read/murray_data/all_gvcf/Euro_gvcfs > /scratch/rjp5nc/Euro_gvcfs.txt
#/project/berglandlab/daphnia_genus/short_read/murray_data/all_gvcf/Euro_gvcfs


cd /scratch/rjp5nc/UK2022_2024/allshortreads/01.RawData/

#need to remove "" from file first
#cut -d',' -f2 "/scratch/rjp5nc/downloadsra.csv" > sra_ids.txt

while read SRA; do
    echo "Processing $SRA..."
    
    # Create a folder for each SRA file
    mkdir -p $SRA

    # Download the SRA file
    prefetch $SRA

    # Convert SRA to FASTQ and compress manually
    fasterq-dump --split-files --outdir "$SRA" "$SRA" && gzip "$SRA"/*.fastq

    echo "$SRA completed!"
done < sra_ids.txt






SRA_ACCESSION="ASM2113471v1"


# Create output directory
OUTPUT_DIR="/scratch/rjp5nc/Reference_genomes/orig_ref"
mkdir -p "$OUTPUT_DIR"

# Download SRA file
echo "Downloading SRA file for $SRA_ACCESSION..."
prefetch --max-size 100G "$SRA_ACCESSION" --output-directory "$OUTPUT_DIR"

# Convert SRA to FASTQ
echo "Converting to FASTQ..."
fasterq-dump --split-files --threads 10 --progress -O "$OUTPUT_DIR" "$OUTPUT_DIR/$SRA_ACCESSION"

# Compress FASTQ files
echo "Compressing FASTQ files..."
gzip "$OUTPUT_DIR/${SRA_ACCESSION}_1.fastq"
gzip "$OUTPUT_DIR/${SRA_ACCESSION}_2.fastq"



