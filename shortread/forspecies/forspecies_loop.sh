
INPUT_FILES=(/scratch/rjp5nc/UK2022_2024/allshortreads/01.RawData/*)

 chmod -x forspecies.sh

# Loop through each input file and submit a job
for INPUT in "${INPUT_FILES[@]}"; do 

    sbatch --job-name=forspecies_$(basename "$INPUT") \
           forspecies.sh "$INPUT"
done




output_file="reference_counts.csv"
echo "Filename,Reference,Count" > "$output_file"  # Create CSV with header

for file in /scratch/rjp5nc/UK2022_2024/allshortreads/counts/*.txt; do  # Loop through all text files
    awk -v filename="$file" '{print filename "," $1 "," $2}' "$file" >> "$output_file"
done

echo "CSV file created: $output_file"



