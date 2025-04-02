#!/usr/bin/env bash
#
#SBATCH -J trimSRA # A single job name for the array
#SBATCH --ntasks-per-node=10 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 2-10:00 # 10 hours
#SBATCH --mem 100G
#SBATCH -o /scratch/rjp5nc/outputerrors/down.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/outputerrors/down.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab


#mv /scratch/rjp5nc/UK2022_2024/allshortreads/01.RawData/SR* /scratch/rjp5nc/UK2022_2024/allshortreads/01.RawData/SRR/

PARENT_DIR="/scratch/rjp5nc/UK2022_2024/allshortreads/01.RawData/SRR/"  # Change this to your actual path
cd "$PARENT_DIR" || exit

for folder in */; do
    folder_name="${folder%/}"  # Remove trailing slash from folder name
    echo "Entering folder: $folder_name"
    cd "$folder_name" || continue  # Enter folder

 if ls *.fastq.gz 1> /dev/null 2>&1; then
        echo "trimming files in $folder_name..."
trimmomatic PE -threads 10 \
${folder_name}_1.fastq.gz \
${folder_name}_2.fastq.gz \
${folder_name}_1.P.trimm.fastq.gz \
${folder_name}_1.U.trimm.fastq.gz \
${folder_name}_2.P.trimm.fastq.gz \
${folder_name}_2.U.trimm.fastq.gz \
ILLUMINACLIP:/home/rjp5nc/miniconda3/bin/trimmomatic/TrimmomaticAdaptors/CombinedPE-PE.fa:2:30:10:8:true

# Merge overlapping reads
/home/rjp5nc/miniconda3/bin/pear \
-f ${folder_name}_1.P.trimm.fastq.gz \
-r ${folder_name}_2.P.trimm.fastq.gz \
-o ${folder_name} \
-j 10
    else
        echo "Warning: No fastq.gz in $folder_name"
    fi
    cd ..
done

