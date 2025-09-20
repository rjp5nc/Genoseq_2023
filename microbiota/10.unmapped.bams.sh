#!/usr/bin/env bash
#SBATCH -J DownloadMap    # Job name
#SBATCH --ntasks=1        # Single task per job
#SBATCH --cpus-per-task=10 # Number of CPU cores per task
#SBATCH -N 1              # Run on one node
#SBATCH -t 1-20:00        # 10 hours runtime
#SBATCH --mem=100G        # Memory per node
#SBATCH -o /scratch/rjp5nc/err/unmapped.%A_%a.out  # Standard output
#SBATCH -e /scratch/rjp5nc/err/unmapped.%A_%a.err  # Standard error
#SBATCH -p standard       # Partition
#SBATCH --account=berglandlab
#SBATCH --array=1-37

module load bwa
module load samtools

# Paths
OUT=/scratch/rjp5nc/UK2022_2024/unmapped_fastqs_newseq
mkdir -p "$OUT/bams" "$OUT/fastqs"
REF=/scratch/rjp5nc/Reference_genomes/post_kraken/US_obtusa_onlydaps.fa

# Get the sample folder for this task
SAMPLE_DIR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" /scratch/rjp5nc/UK2022_2024/unmapped_fastqs_newseq/epistylis.txt)
sample=$(basename "$SAMPLE_DIR")

# Find paired FASTQs inside the folder
R1=$(ls /scratch/rjp5nc/UK2022_2024/allshortreads/01.RawData/newseq/"$SAMPLE_DIR"/*_trimmedmerged1.fq.gz)
R2=$(ls /scratch/rjp5nc/UK2022_2024/allshortreads/01.RawData/newseq/"$SAMPLE_DIR"/*_trimmedmerged2.fq.gz)

echo "[$SLURM_ARRAY_TASK_ID] Mapping $sample ..."

# Map and produce BAM
bwa mem -t 10 "$REF" "$R1" "$R2" \
    | samtools view -bS -o "$OUT/bams/${sample}.bam"

# Extract unmapped reads
samtools fastq -f 4 "$OUT/bams/${sample}.bam" \
    -1 "$OUT/fastqs/${sample}_R1.unmapped.fq.gz" \
    -2 "$OUT/fastqs/${sample}_R2.unmapped.fq.gz" \
    -0 /dev/null -s /dev/null -n

echo "[$SLURM_ARRAY_TASK_ID] Done $sample"