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
PARENT_DIR="/scratch/rjp5nc/UK2022_2024/allshortreads/01.RawData/newseq/"  # Change this to your actual path
cd "$PARENT_DIR" || exit

# Loop through all subdirectories
for folder in */; do
    folder_name="${folder%/}"  # Remove trailing slash from folder name
    echo "Entering folder: $folder_name"
    cd "$folder_name" || continue  # Enter folder

    # Check for forward reads
    if ls *trimmed1.fq.gz 1> /dev/null 2>&1; then
        echo "Merging all forward reads in $folder_name..."
        cat *trimmed1.fq.gz > "${folder_name}_trimmedmerged1.fq.gz"
        rm *trimmed1.fq.gz
        rm *_1.fq.gz
    else
        echo "Warning: No fq.gz in $folder_name"
    fi

    # Check for reverse reads
    if ls *trimmed2.fq.gz 1> /dev/null 2>&1; then
        echo "Merging all reverse reads in $folder_name..."
        cat *trimmed2.fq.gz > "${folder_name}_trimmedmerged2.fq.gz"
        rm *trimmed2.fq.gz
        rm *_2.fq.gz
    else
        echo "Warning: No reverse reads found in $folder_name"
    fi

    cd ..  # Move back to parent directory
done
echo "Done merging"












#for folder in */; do
 #   folder_name="${folder%/}"  # Remove trailing slash from folder name
  #  echo "Entering folder: $folder_name"
   # cd "$folder_name" || continue  # Enter folder

    # Check for forward reads
    #if ls *.fastq.gz 1> /dev/null 2>&1; then
     #   echo "Merging all forward reads in $folder_name..."
# Trim illumina adaptors



