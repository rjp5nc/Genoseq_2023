#!/usr/bin/env bash
#
#SBATCH -J trimSRA # A single job name for the array
#SBATCH --ntasks-per-node=10 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-10:00 # 10 hours
#SBATCH --mem 100G
#SBATCH -o /scratch/rjp5nc/erroroutputs/down.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/erroroutputs/down.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

###sbatch --array=1-287 510.trim_SRR_fastq_mito.sh

PARENT_DIR="/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/fasta"

SRA_LIST="/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/sra_USobtusa.txt"
BASE_OUTDIR="/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/fasta"

cd /scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/fasta

# Get SRA accession for this array task
SRA=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SRA_LIST")
echo $SRA
R1="${SRA}_1.fastq.gz"
R2="${SRA}_2.fastq.gz"


    # Run Trimmomatic
    trimmomatic PE -threads 10 \
        ${SRA}_1.fastq.gz \
        ${SRA}_2.fastq.gz \
        ${SRA}_1.P.trimm.fastq.gz \
        ${SRA}_1.U.trimm.fastq.gz \
        ${SRA}_2.P.trimm.fastq.gz \
        ${SRA}_2.U.trimm.fastq.gz \
        ILLUMINACLIP:/home/rjp5nc/miniconda3/bin/trimmomatic/TrimmomaticAdaptors/CombinedPE-PE.fa:2:30:10:8:true

    # Run PEAR to merge overlapping reads
    /home/rjp5nc/miniconda3/bin/pear \
        -f ${SRA}_1.P.trimm.fastq.gz \
        -r ${SRA}_2.P.trimm.fastq.gz \
        -o ${SRA} \
        -j 10

mv /scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/fasta/${SRA}.assembled.fastq /scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/fasta_assembled
mv /scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/fasta/${SRA}.unassembled*.fastq /scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/fasta_assembled
mv /scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/fasta/${SRA}.disc*.fastq /scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/fasta_assembled

rm /scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/fasta/${SRA}*

