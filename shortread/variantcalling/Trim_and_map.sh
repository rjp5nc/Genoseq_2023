#!/usr/bin/env bash
#
#SBATCH -J download # A single job name for the array
#SBATCH --ntasks-per-node=10 # ten cores
#SBATCH -N 1 # on one node
#SBATCH -t 16:00:00 ### 8 hours
#SBATCH --mem 60G
#SBATCH -o /scratch/rjp5nc/outputerrors/download_%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/outputerrors/download_%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

module load cutadapt gcc/11.4.0 bwa/0.7.17 samtools/1.17 picard/2.27.5

# sbatch --array=2 ~/Genoseq_2023/Trim_and_map.sh
### sacct -j 40641767
###  cat /scratch/rjp5nc/download_40446865_1.err
# sacct -j

### example:
# SLURM_ARRAY_TASK_ID=1

###dir=/project/berglandlab/Robert/UKSequencing2022_2024/usftp21.novogene.com/01.RawData/Ro*

dir=$1
echo ${dir}

L6_1=$( ls ${dir}/*L6_1.fq.gz )
L6_2=$( ls ${dir}/*L6_2.fq.gz )
L3_1=$( ls ${dir}/*L3_1.fq.gz )
L3_2=$( ls ${dir}/*L3_2.fq.gz )
L4_1=$( ls ${dir}/*L4_1.fq.gz )
L4_2=$( ls ${dir}/*L4_2.fq.gz )


sample_name=$(basename ${dir})


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
${L6_1}.trimmed1.fq.gz ${L6_2}.trimmed2.fq.gz |
samtools view -@ 10 -Sbh -q 20 -F 0x100 - > ${dir}/${sample_name}.L6.bam


cutadapt \
-q 18 \
--minimum-length 75 \
-o ${L3_1}.trimmed1.fq.gz \
-p ${L3_2}.trimmed2.fq.gz \
-O 15 \
-n 3 \
--cores=10 \
${L3_1} ${L3_2}

bwa mem \
-t 10 \
-R "@RG\tID:${dir}_L3\tSM:sample_name\tPL:illumina\tLB:lib1" \
/project/berglandlab/daphnia_ref/totalHiCwithallbestgapclosed.fa \
${L3_1}.trimmed1.fq.gz ${L3_2}.trimmed2.fq.gz |
samtools view -@ 10 -Sbh -q 20 -F 0x100 - > ${dir}/${sample_name}.L3.bam

### run for L4
cutadapt \
-q 18 \
--minimum-length 75 \
-o ${L4_1}.trimmed1.fq.gz \
-p ${L4_2}.trimmed2.fq.gz \
-O 15 \
-n 3 \
--cores=10 \
${L4_1} ${L4_2}

bwa mem \
-t 10 \
-R "@RG\tID:${dir}_L4\tSM:sample_name\tPL:illumina\tLB:lib1" \
/project/berglandlab/daphnia_ref/totalHiCwithallbestgapclosed.fa \
${L4_1}.trimmed1.fq.gz ${L4_2}.trimmed2.fq.gz |
samtools view -@ 10 -Sbh -q 20 -F 0x100 - > ${dir}/${sample_name}.L4.bam


java -jar $EBROOTPICARD/picard.jar MergeSamFiles \
    -I ${dir}/${sample_name}.L3.bam \
    -I ${dir}/${sample_name}.L4.bam \
    -I ${dir}/${sample_name}.L6.bam \
    -O /scratch/rjp5nc/UK2022_2024/allshortreads/mergedbams/${sample_name}.merged.bam \
    --SORT_ORDER unsorted


java -jar $EBROOTPICARD/picard.jar SortSam \
    -I /scratch/rjp5nc/UK2022_2024/allshortreads/mergedbams/${sample_name}.merged.bam \
    -O /scratch/rjp5nc/UK2022_2024/allshortreads/sortedbams/${sample_name}.sorted.bam \
    -SORT_ORDER coordinate

### PCR duplicate removal
java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
-REMOVE_DUPLICATES true \
-I /scratch/rjp5nc/UK2022_2024/allshortreads/sortedbams/${sample_name}.sorted.bam \
-O /scratch/rjp5nc/UK2022_2024/allshortreads/sortedbamsdedup/${sample_name}.sort.dedup.bam \
-M /scratch/rjp5nc/UK2022_2024/allshortreads/sortedbamsreport/${sample_name}.mark_duplicates_report.txt \
-VALIDATION_STRINGENCY SILENT