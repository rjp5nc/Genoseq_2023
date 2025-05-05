#!/usr/bin/env bash
#
#SBATCH -J forspecies # A single job name for the array
#SBATCH --ntasks-per-node=10 # ten core
#SBATCH -N 1 # on one node
#SBATCH -t 0-10:00  ### 10 hours
#SBATCH --mem 60G
#SBATCH -o /scratch/rjp5nc/error_out/forspecies.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/error_out/forspecies.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab
 
module load cutadapt gcc/11.4.0 bwa/0.7.17 samtools/1.17 picard/2.27.5

dir=$1
echo ${dir}

sample_name=$(basename ${dir})

L6_1=$( ls ${dir}/*L6_1.fq.gz )
L6_2=$( ls ${dir}/*L6_2.fq.gz )
L3_1=$( ls ${dir}/*L3_1.fq.gz )
L3_2=$( ls ${dir}/*L3_2.fq.gz )
L4_1=$( ls ${dir}/*L4_1.fq.gz )
L4_2=$( ls ${dir}/*L4_2.fq.gz )

#Run once
#bwa index /scratch/rjp5nc/UK2022_2024/allshortreads/allcoi_updated22.fa


outputdir="/scratch/rjp5nc/UK2022_2024/allshortreads/outputcoibams"

bwa mem \
-t 10 \
-R "@RG\tID:${dir}_L6\tSM:sample_name\tPL:illumina\tLB:lib1" \
/scratch/rjp5nc/UK2022_2024/allshortreads/allcoi_updated22.fa \
${L6_1}.trimmed1.fq.gz ${L6_2}.trimmed2.fq.gz |
samtools view -@ 10 -Sbh -q 20 -F 0x100 - > ${outputdir}/${sample_name}.L6.bam

#samtools view ${outputdir}/${sample_name}.L6.bam | less -S

bwa mem \
-t 10 \
-R "@RG\tID:${dir}_L3\tSM:sample_name\tPL:illumina\tLB:lib1" \
/scratch/rjp5nc/UK2022_2024/allshortreads/allcoi_updated22.fa \
${L3_1}.trimmed1.fq.gz ${L3_2}.trimmed2.fq.gz |
samtools view -@ 10 -Sbh -q 20 -F 0x100 - > ${outputdir}/${sample_name}.L3.bam

bwa mem \
-t 10 \
-R "@RG\tID:${dir}_L4\tSM:sample_name\tPL:illumina\tLB:lib1" \
/scratch/rjp5nc/UK2022_2024/allshortreads/allcoi_updated22.fa \
${L4_1}.trimmed1.fq.gz ${L4_2}.trimmed2.fq.gz |
samtools view -@ 10 -Sbh -q 20 -F 0x100 - > ${outputdir}/${sample_name}.L4.bam


java -jar $EBROOTPICARD/picard.jar MergeSamFiles \
    -I ${outputdir}/${sample_name}.L3.bam \
    -I ${outputdir}/${sample_name}.L4.bam \
    -I ${outputdir}/${sample_name}.L6.bam \
    -O /scratch/rjp5nc/UK2022_2024/allshortreads/mitobams/${sample_name}.merged.bam \
    --SORT_ORDER unsorted

java -jar $EBROOTPICARD/picard.jar SortSam \
    -I /scratch/rjp5nc/UK2022_2024/allshortreads/mitobams/${sample_name}.merged.bam \
    -O /scratch/rjp5nc/UK2022_2024/allshortreads/mitosortedbams/${sample_name}.sorted.bam \
    -SORT_ORDER coordinate

samtools view /scratch/rjp5nc/UK2022_2024/allshortreads/P2_2.bam | awk '{count[$3]++} END {for (val in count) print val, count[val]}' | sort -k2,2nr >> /scratch/rjp5nc/UK2022_2024/allshortreads/P2_2counts.txt
