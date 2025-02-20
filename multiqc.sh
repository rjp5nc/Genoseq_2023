#!/usr/bin/env bash
#
#SBATCH -J Multiqc # A single job name for the array
#SBATCH --ntasks-per-node=40 # ten cores
#SBATCH -N 1 # on one node
#SBATCH -t 24:00:00 ### 8 hours
#SBATCH --mem 200G
#SBATCH -o /scratch/rjp5nc/Multiqc.out # Standard output
#SBATCH -e /scratch/rjp5nc/Multiqc.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab


# sbatch ~/Genoseq_2023/multiqc.sh
### sacct -j 45333345
###  cat /scratch/rjp5nc/download_40446865_1.err
# sacct -j  
module load fastqc
module load multiqc

INPUT_DIR="/scratch/rjp5nc/UK2022_2024/allshortreads/01.RawData/"
OUTPUT_DIR="/scratch/rjp5nc/outputdir2lanes"

find "$INPUT_DIR" -type f -name "*L3_1.fq.gz" | while read -r FILE; do
    echo "Running FastQC on $FILE..."
    fastqc -o "$OUTPUT_DIR" "$FILE" --no-graphical
done

#multiqc -o /scratch/rjp5nc/multiqc2lanes /scratch/rjp5nc/outputdir2lanes