#!/usr/bin/env bash
#
#SBATCH -J download # A single job name for the array
#SBATCH --ntasks-per-node=10 # ten cores
#SBATCH -N 1 # on one node
#SBATCH -t 16:00:00 ### 8 hours
#SBATCH --mem 60G
#SBATCH -o /scratch/rjp5nc/download_%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/download_%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

module load cutadapt gcc/11.4.0 bwa/0.7.17 samtools/1.17 picard/2.27.5

# sbatch --array=2 ~/Genoseq_2023/Trim_and_map.sh
### sacct -j 40641767
###  cat /scratch/rjp5nc/download_40446865_1.err
# sacct -j

### example:
# SLURM_ARRAY_TASK_ID=1

dir=/project/berglandlab/Robert/UKSequencing2022_2024/usftp21.novogene.com/01.RawData/Gilmer5_A1





java -jar $EBROOTPICARD/picard.jar SortSam \
-I ${dir}.L6.bam \
-O ${dir}.L6.sorted.bam \
-SORT_ORDER coordinate

### PCR duplicate removal
java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
-REMOVE_DUPLICATES true \
-I ${dir}.L6.sorted.bam \
-O ${dir}.sort.dedup.bam \
-M ${dir}.mark_duplicates_report.txt \
-VALIDATION_STRINGENCY SILENT