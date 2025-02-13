
#!/usr/bin/env bash
#
#SBATCH -J trimloop # A single job name for the array
#SBATCH --ntasks-per-node=10 # ten cores
#SBATCH -N 1 # on one node
#SBATCH -t 16:00:00 ### 8 hours
#SBATCH --mem 60G
#SBATCH -o /scratch/rjp5nc/outputerrors/download_%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/outputerrors/download_%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

INPUT_FILES=(/scratch/rjp5nc/UK2022_2024/allshortreads/01.RawData/*)

 chmod -x Trim_and_map.sh

# Loop through each input file and submit a job
for INPUT in "${INPUT_FILES[@]}"; do 

    sbatch --job-name=TrimMap_$(basename "$INPUT") \
           Trim_and_map.sh "$INPUT"
done


#For indexing

