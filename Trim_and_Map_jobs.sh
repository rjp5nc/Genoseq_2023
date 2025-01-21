


INPUT_FILES=("/project/berglandlab/Robert/UKSequencing2022_2024/usftp21.novogene.com/01.RawData/Gilmer5_A1")  # List of input file prefixes (without extensions)

chmod -x Trim_and_map.sh

# Loop through each input file and submit a job
for INPUT in "${INPUT_FILES[@]}"; do 
    sbatch --job-name=TrimMap_${INPUT} \
           Trim_and_map.sh ${INPUT}
done

sacct -j 1466105






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
samtools view -@ 10 -Sbh -q 20 -F 0x100 - > ${dir}.L6.bam
### samtools view -@ 10 -Sbh -q 20 -F 0x100 - > ${dir}/${dir}.L6.bam
