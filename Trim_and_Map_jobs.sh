


INPUT_FILES=(/project/berglandlab/Robert/UKSequencing2022_2024/usftp21.novogene.com/01.RawData/Rob*)

### chmod -x Trim_and_map.sh

# Loop through each input file and submit a job
for INPUT in "${INPUT_FILES[@]}"; do 

    sbatch --job-name=TrimMap_$(basename "$INPUT") \
           Trim_and_map.sh "$INPUT"
done

sacct -j 1466105


for INPUT in "${INPUT_FILES[@]}"; do 
    sbatch --job-name=TrimMap_${INPUT} \
           Trim_and_map.sh ${INPUT}
done


