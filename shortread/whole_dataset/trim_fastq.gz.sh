#!/usr/bin/env bash
#
#SBATCH -J trimSRA # A single job name for the array
#SBATCH --ntasks-per-node=10 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-10:00 # 10 hours
#SBATCH --mem 100G
#SBATCH -o /scratch/rjp5nc/outputerrors/down.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/outputerrors/down.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab
#SBATCH --array=1-576%20   # Adjust the range based on the number of folders

#mv /scratch/rjp5nc/UK2022_2024/allshortreads/01.RawData/SR* /scratch/rjp5nc/UK2022_2024/allshortreads/01.RawData/SRR/

PARENT_DIR="/scratch/rjp5nc/UK2022_2024/allshortreads/01.RawData/newseq"

# Get list of all folders
FOLDERS=($(ls -d ${PARENT_DIR}/*/))  # Array of folder paths

# Select the folder based on the job array index
FOLDER="${FOLDERS[$SLURM_ARRAY_TASK_ID - 1]}"

# Extract folder name
FOLDER_NAME=$(basename "$FOLDER")

# Change to the working directory
cd "$FOLDER" || exit 1

# Check if fastq.gz files exist
if ls *.fq.gz 1> /dev/null 2>&1; then
    echo "Processing $FOLDER_NAME"

    # Run Trimmomatic
 #   trimmomatic PE -threads 10 \
 #       ${FOLDER_NAME}_trimmedmerged1.fq.gz \
 #       ${FOLDER_NAME}_trimmedmerged1.fq.gz \
 #       ${FOLDER_NAME}_1.P.trimm.fastq.gz \
 #       ${FOLDER_NAME}_1.U.trimm.fastq.gz \
 #       ${FOLDER_NAME}_2.P.trimm.fastq.gz \
 #       ${FOLDER_NAME}_2.U.trimm.fastq.gz \
 #       ILLUMINACLIP:/home/rjp5nc/miniconda3/bin/trimmomatic/TrimmomaticAdaptors/CombinedPE-PE.fa:2:30:10:8:true

    # Run PEAR to merge overlapping reads
    /home/rjp5nc/miniconda3/bin/pear \
        -f ${FOLDER_NAME}_trimmedmerged1.fq.gz \
        -r ${FOLDER_NAME}_trimmedmerged2.fq.gz \
        -o ${FOLDER_NAME} \
        -j 10
else
    echo "Warning: No fq.gz files in $FOLDER_NAME"
fi