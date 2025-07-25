#!/usr/bin/env bash
#SBATCH -J DownloadMap    # Job name
#SBATCH --ntasks=1        # Single task per job
#SBATCH --cpus-per-task=10 # Number of CPU cores per task
#SBATCH -N 1              # Run on one node
#SBATCH -t 1-20:00        # 10 hours runtime
#SBATCH --mem=100G        # Memory per node
#SBATCH -o /scratch/rjp5nc/outputerrors/down.%A_%a.out  # Standard output
#SBATCH -e /scratch/rjp5nc/outputerrors/down.%A_%a.err  # Standard error
#SBATCH -p standard       # Partition
#SBATCH --account=berglandlab
#SBATCH --array=1-21   # Adjust based on the number of samples

# Load necessary modules
module load gcc htslib
module load sratoolkit/3.1.1
module load trimmomatic
module load bwa
module load samtools
module load picard

#SRR5012770
#SRR5012773
# 3387565_19 

# Define working directories
infq="/scratch/rjp5nc/UK2022_2024/allshortreads/01.RawData/SRR"
outfq="/scratch/rjp5nc/UK2022_2024/allshortreads/01.RawData/SRR_mito2"
outbam="/scratch/rjp5nc/UK2022_2024/mapped_bam_newmito_newref"

# Ensure output directories exist
mkdir -p "${outfq}" "${outbam}" "${infq}"

# Read sample ID and reference path from CSV
CSV_FILE="/scratch/rjp5nc/UK2022_2024/forref2_mito.csv"
line=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" ${CSV_FILE})

# Extract fields (assuming CSV format: sample_id,reference_path)
samp=$(echo "$line" | cut -d',' -f4)
ref_path=$(echo "$line" | cut -d',' -f5)
ref_path=$(echo "${ref_path}" | tr -d '\r')

#samp=SRR14370492
#ref_path=/scratch/rjp5nc/Reference_genomes/post_kraken/assembly.hap2_onlydaps.fasta

# Ensure sample and reference path are valid
if [[ -z "$samp" || -z "$ref_path" ]]; then
    echo "Error: Invalid sample or reference genome for task ID ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

echo "Processing sample: ${samp} with reference: ${ref_path}"

chmod u+w /scratch/rjp5nc/UK2022_2024/allshortreads/01.RawData/SRR/${samp}/${samp}.*.fastq

# Map to reference genome (assembled reads)
bwa mem -t 10 -K 100000000 -Y ${ref_path} ${infq}/${samp}/${samp}.assembled.fastq | \
samtools view -Suh -q 20 -F 0x100 | \
samtools sort --threads 10 -o ${outfq}/${samp}.sort.bam
samtools index ${outfq}/${samp}.sort.bam

# Map unassembled reads
bwa mem -t 10 -K 100000000 -Y ${ref_path} \
${infq}/${samp}/${samp}.unassembled.forward.fastq  \
${infq}/${samp}/${samp}.unassembled.reverse.fastq | \
samtools view -Suh -q 20 -F 0x100 | \
samtools sort --threads 10 -o ${outfq}/${samp}.filt.unassembled.sort.bam
samtools index ${outfq}/${samp}.filt.unassembled.sort.bam

# Merge assembled and unassembled BAM files
samtools merge ${outfq}/${samp}.filt.merged.bam \
    ${outfq}/${samp}.sort.bam \
    ${outfq}/${samp}.filt.unassembled.sort.bam

# Index merged BAM
samtools index ${outfq}/${samp}.filt.merged.bam

# Mark duplicates
java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
    REMOVE_DUPLICATES=true \
    INPUT=${outfq}/${samp}.filt.merged.bam \
    OUTPUT=${outfq}/${samp}_finalmap.bam \
    METRICS_FILE=${outfq}/${samp}_finalmap_mdups.metrics \
    CREATE_INDEX=true

# Move final BAM to output directory
cp ${outfq}/${samp}_finalmap* ${outbam}/

# Remove intermediate files
#rm -f ${outfq}/${samp}.*

echo "Finished processing ${samp}"