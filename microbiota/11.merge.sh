#!/usr/bin/env bash
#
#SBATCH -J merge # A single job name for the array
#SBATCH --ntasks-per-node=10 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 2-10:00 # 10 hours
#SBATCH --mem 100G
#SBATCH -o /scratch/rjp5nc/erroroutputs/down.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/erroroutputs/down.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

#cp -r /scratch/rjp5nc/UK2022_2024/allshortreads/01.RawData/Rockpool2_B3/ /scratch/rjp5nc/UK2022_2024/allshortreads/01.RawData/testing/
# Set the parent directory containing all the folders
PARENT_DIR="/scratch/rjp5nc/UK2022_2024/unmapped_fastqs_newseq/fastqs"  # Change this to your actual path
cd "$PARENT_DIR" || exit

# Loop through all subdirectories
# Check for forward reads
if ls *R1.unmapped.fq.gz 1> /dev/null 2>&1; then
    cat *R1.unmapped.fq.gz > "unmapped_trimmedmerged1.fq.gz"
fi

# Check for reverse reads
if ls *R2.unmapped.fq.gz 1> /dev/null 2>&1; then
    cat *R2.unmapped.fq.gz > "unmapped_trimmedmerged2.fq.gz"
fi








#for folder in */; do
 #   folder_name="${folder%/}"  # Remove trailing slash from folder name
  #  echo "Entering folder: $folder_name"
   # cd "$folder_name" || continue  # Enter folder

    # Check for forward reads
    #if ls *.fastq.gz 1> /dev/null 2>&1; then
     #   echo "Merging all forward reads in $folder_name..."
# Trim illumina adaptors



