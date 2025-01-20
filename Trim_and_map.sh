#!/usr/bin/env bash
#
#SBATCH -J download # A single job name for the array
#SBATCH --ntasks-per-node=10 # ten cores
#SBATCH -N 1 # on one node
#SBATCH -t 8:00:00 ### 8 hours
#SBATCH --mem 60G
#SBATCH -o /scratch/rjp5nc/download_%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/download_%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

module load cutadapt gcc/11.4.0 bwa/0.7.17 samtools/1.17
java -jar $EBROOTPICARD/picard.jar

# sbatch --array=2 ~/Genoseq_2023/Trim_and_map.sh
### sacct -j 40641767
###  cat /scratch/rjp5nc/download_40446865_1.err
# sacct -j

### example:
# SLURM_ARRAY_TASK_ID=1

dir=/project/berglandlab/Robert/UKSequencing2022_2024/usftp21.novogene.com/01.RawData/Gilmer5_A1
L6_1=$( ls ${dir}/*L6_1.fq.gz )
L6_2=$( ls ${dir}/*L6_2.fq.gz )

## adapter removal

#### adapter removal and inital mapping
### run for L3
cutadapt \
-q 18 \
--minimum-length 75 \
-o ${L6_1}.trimmed1.fq.gz \
-p ${L6_2}.trimmed2.fq.gz \
-O 15 \
-n 3 \
--cores=10 \
${L6_1} ${L6_2}

bwa mem \
-t 10 \
-R "@RG\tID:${dir}_L6\tSM:sample_name\tPL:illumina\tLB:lib1" \
/project/berglandlab/daphnia_ref/totalHiCwithallbestgapclosed.fa \
${L6_1} ${L6_2} |
samtools view -@ 10 -Sbh -q 20 -F 0x100 - > /Bams/${dir}.L6.bam
### samtools view -@ 10 -Sbh -q 20 -F 0x100 - > ${dir}/${dir}.L6.bam

### PCR duplicate removal
java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
-REMOVE_DUPLICATES true \
-I /Bams/${dir}.L6.bam \
-O /Bams/${dir}.sort.dedup.bam \
-M /Bams/${dir}.mark_duplicates_report.txt \
-VALIDATION_STRINGENCY SILENT
