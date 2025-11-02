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
/scratch/rjp5nc/UK2022_2024/allshortreads/01.RawData/newseq/RobertUK_G5/RobertUK_G5.assembled.fastq
/scratch/rjp5nc/UK2022_2024/allshortreads/01.RawData/newseq/RobertUK_G5/RobertUK_G5.unassembled.forward.fastq
/scratch/rjp5nc/UK2022_2024/allshortreads/01.RawData/newseq/RobertUK_G5/RobertUK_G5.unassembled.reverse.fastq

samp=Canal1Reed
# Map to reference genome (assembled reads)
bwa mem -t 10 -K 100000000 -Y ${ref_path} /scratch/rjp5nc/UK2022_2024/allshortreads/01.RawData/newseq/RobertUK_G12/RobertUK_G12.assembled.fastq | \
samtools view -Suh -q 20 -F 0x100 | \
samtools sort --threads 10 -o ${outfq}/chlorella_Canal1Reed.sort.bam
samtools index ${outfq}/chlorella_Canal1Reed.sort.bam

# Map unassembled reads
bwa mem -t 10 -K 100000000 -Y ${ref_path} \
/scratch/rjp5nc/UK2022_2024/allshortreads/01.RawData/newseq/RobertUK_G12/RobertUK_G12.unassembled.forward.fastq  \
/scratch/rjp5nc/UK2022_2024/allshortreads/01.RawData/newseq/RobertUK_G12/RobertUK_G12.unassembled.reverse.fastq | \
samtools view -Suh -q 20 -F 0x100 | \
samtools sort --threads 10 -o ${outfq}/${samp}.filt.unassembled.sort.bam
samtools index ${outfq}/${samp}.filt.unassembled.sort.bam

# Merge assembled and unassembled BAM files
samtools merge ${outfq}/${samp}.filt.merged.bam \
    ${outfq}/chlorella_Canal1Reed.sort.bam \
    ${outfq}/${samp}.filt.unassembled.sort.bam

# Index merged BAM
samtools index ${outfq}/${samp}.filt.merged.bam

# Mark duplicates
java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
    REMOVE_DUPLICATES=true \
    INPUT=${outfq}/${samp}.filt.merged.bam \
    OUTPUT=${outfq}/${samp}.filt.merged.bam \
    METRICS_FILE=${outfq}/${samp}_finalmap_mdups.metrics \
    CREATE_INDEX=true

# Move final BAM to output directory
mv ${outfq}/${samp}.filt.merged.* ${outbam}/

# Remove intermediate files
rm -f ${outfq}/${samp}.*

echo "Finished processing ${samp}"

#RobertUKG5 for Reed


module load gatk


samtools faidx /scratch/rjp5nc/Reference_genomes/chlorella_ref/GCA_023343905.1_cvul_genomic.fa

gatk CreateSequenceDictionary \
   -R /scratch/rjp5nc/Reference_genomes/chlorella_ref/GCA_023343905.1_cvul_genomic.fa \
   -O /scratch/rjp5nc/Reference_genomes/chlorella_ref/GCA_023343905.1_cvul_genomic.dict

java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
-I /scratch/rjp5nc/microbiota/chlorella/fastq/${samp}.filt.merged.bam \
-O /scratch/rjp5nc/microbiota/chlorella/fastq/${samp}.finalmap_rg.bam \
-LB "library" \
-PL "ILLumina" \
-PU "platunit" \
-SM chlorella_reed_canal1

# Index Bam files
java -jar $EBROOTPICARD/picard.jar BuildBamIndex \
-I /scratch/rjp5nc/microbiota/chlorella/fastq/${samp}.finalmap_rg.bam



#Need to make gvcf files of more than just one bam file

gatk HaplotypeCaller \
-R ${ref_path} \
-I /scratch/rjp5nc/microbiota/chlorella/fastq/${samp}.finalmap_rg.bam \
-O /scratch/rjp5nc/microbiota/chlorella/gvcf/chlorella_Reed_canal1.g.vcf \
-ERC GVCF

gatk HaplotypeCaller \
-R ${ref_path} \
-I /scratch/rjp5nc/microbiota/chlorella/mapped_bam/cholerlla_ours_finalmap_RG.bam \
-O /scratch/rjp5nc/microbiota/chlorella/gvcf/chlorella_ours.g.vcf \
-ERC GVCF



gatk CombineGVCFs \
   -R ${ref_path} \
   --variant /scratch/rjp5nc/microbiota/chlorella/gvcf/chlorella_ours.g.vcf  \
   --variant  /scratch/rjp5nc/microbiota/chlorella/gvcf/chlorella_Reed_canal1.g.vcf  \
   -O  /scratch/rjp5nc/microbiota/chlorella/gvcf/chlorella.g.vcf 

module load gcc/14.2.0 htslib/1.17

# Compress and index with Tabix
bgzip /scratch/rjp5nc/microbiota/chlorella/gvcf/chlorella.g.vcf 
tabix -p vcf /scratch/rjp5nc/microbiota/chlorella/gvcf/chlorella.g.vcf.gz

# Finish



JAVAMEM=40g
gatk --java-options "-Xmx${JAVAMEM}" GenotypeGVCFs \
   -R ${ref_path} \
   -V /scratch/rjp5nc/microbiota/chlorella/gvcf/chlorella.g.vcf.gz \
   -O /scratch/rjp5nc/microbiota/chlorella/gvcf/chlorella_genotype.vcf.gz
