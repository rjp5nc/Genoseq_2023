#!/usr/bin/env bash
#SBATCH -J makebams    # Job name
#SBATCH --ntasks=1        # Single task per job
#SBATCH --cpus-per-task=10 # Number of CPU cores per task
#SBATCH -N 1              # Run on one node
#SBATCH -t 0-10:00        # 10 hours runtime
#SBATCH --mem=100G        # Memory per node
#SBATCH -o /scratch/rjp5nc/erroroutputs/down.%A_%a.out  # Standard output
#SBATCH -e /scratch/rjp5nc/erroroutputs/down.%A_%a.err  # Standard error
#SBATCH -p standard       # Partition
#SBATCH --account=berglandlab
#SBATCH --array=1-576%40   # Adjust based on the number of samples

# Load necessary modules
module load gcc htslib
module load sratoolkit/3.1.1
module load trimmomatic
module load bwa
module load samtools
module load picard

# Define working directories
outfq="/scratch/rjp5nc/microbiota/chlorella/fastq"
outbam="/scratch/rjp5nc/microbiota/chlorella/mapped_bam"

# Ensure output directories exist
mkdir -p "${outfq}" "${outbam}"

# Extract fields (assuming CSV format: sample_id,reference_path)
ref_path=/scratch/rjp5nc/Reference_genomes/chlorella_ref/GCA_023343905.1_cvul_genomic.fa

#ref_path=/scratch/rjp5nc/Reference_genomes/post_kraken/assembly.hap2_onlydaps.fasta

# Map to reference genome (assembled reads)
bwa mem -t 10 -K 100000000 -Y ${ref_path} /scratch/rjp5nc/HMW/HMWDNAElvis3/m84128_250121_222443_s2.hifi_reads.bc2104.fastq | \
samtools view -Suh -q 20 -F 0x100 | \
samtools sort --threads 10 -o ${outfq}/chlorella_Reed.sort.bam
samtools index ${outfq}/${samp}.sort.bam

# Map unassembled reads
bwa mem -t 10 -K 100000000 -Y ${ref_path} \
${outfq}/${samp}/${samp}.unassembled.forward.fastq  \
${outfq}/${samp}/${samp}.unassembled.reverse.fastq | \
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
    INPUT=${outfq}/chlorella_Reed.sort.bam \
    OUTPUT=${outfq}/chlorella_Reed_finalmap.bam \
    METRICS_FILE=${outfq}/${samp}_finalmap_mdups.metrics \
    CREATE_INDEX=true

# Move final BAM to output directory
mv ${outfq}/${samp}_finalmap* ${outbam}/

# Remove intermediate files
rm -f ${outfq}/${samp}.*

echo "Finished processing ${samp}"




module load gatk


samtools faidx /scratch/rjp5nc/Reference_genomes/chlorella_ref/GCA_023343905.1_cvul_genomic.fa

gatk CreateSequenceDictionary \
   -R /scratch/rjp5nc/Reference_genomes/chlorella_ref/GCA_023343905.1_cvul_genomic.fa \
   -O /scratch/rjp5nc/Reference_genomes/chlorella_ref/GCA_023343905.1_cvul_genomic.dict

java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
-I /scratch/rjp5nc/microbiota/chlorella/fastq/chlorella_Reed_finalmap.bam \
-O /scratch/rjp5nc/microbiota/chlorella/fastq/chlorella_Reed_finalmap_RG.bam \
-LB "library" \
-PL "ILLumina" \
-PU "platunit" \
-SM chlorella_reed

# Index Bam files
java -jar $EBROOTPICARD/picard.jar BuildBamIndex \
-I /scratch/rjp5nc/microbiota/chlorella/fastq/chlorella_Reed_finalmap_RG.bam



#Need to make gvcf files of more than just one bam file

gatk HaplotypeCaller \
-R ${ref_path} \
-I /scratch/rjp5nc/microbiota/chlorella/fastq/chlorella_Reed_finalmap_RG.bam \
-O /scratch/rjp5nc/microbiota/chlorella/fastq/chlorella_Reed.g.vcf \
-ERC GVCF


