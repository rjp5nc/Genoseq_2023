
INPUT_FILES=(/scratch/rjp5nc/UK2022_2024/allshortreads/01.RawData/*)

 chmod -x forspecies.sh

# Loop through each input file and submit a job
for INPUT in "${INPUT_FILES[@]}"; do 

    sbatch --job-name=forspecies_$(basename "$INPUT") \
           Trim_and_map.sh "$INPUT"
done
