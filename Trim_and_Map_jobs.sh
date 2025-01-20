


INPUT_FILES=("/project/berglandlab/Robert/UKSequencing2022_2024/usftp21.novogene.com/01.RawData/Gilmer5_A1")  # List of input file prefixes (without extensions)

chmod -x Trim_and_Map.sh

# Loop through each input file and submit a job
for INPUT in "${INPUT_FILES[@]}"; do 
    sbatch --job-name=TrimMap_${INPUT} \
           Trim_and_map.sh ${INPUT}
done