#!/usr/bin/env bash
#
#SBATCH -J windowfasta # A single job name for the array
#SBATCH --ntasks-per-node=5 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-2:00:00 # 8 hours
#SBATCH --mem 60G
#SBATCH -o /scratch/rjp5nc/err/gatk.chrom.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/err/gatk.chrom.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab


# Input FASTA file
FASTA_FILE="/scratch/rjp5nc/Reference_genomes/post_kraken/assembly.hap2_onlydaps.fasta"
OUTPUT_FILE="/scratch/rjp5nc/UK2022_2024/daphnia_phylo/interval_DBI_paramList_euobtusa.txt"

# Set window size (1000000)
WINDOW_SIZE=1000000

# Write header to output file
echo "id,seqnames,start,end,width" > $OUTPUT_FILE

# Initialize variables
id=1

# Read the FASTA file line by line
while IFS= read -r line; do
    # Check if the line is a FASTA header (starts with '>')
    if [[ ${line} == ">"* ]]; then
        # Extract the sequence name from the header (remove the '>' character)
        seqname=${line:1}
        seq=""
    else
        # Append the sequence to a variable (the entire sequence for this scaffold)
        seq+=$line
    fi
done < "$FASTA_FILE"

# Now process the sequence for the scaffold
len=${#seq}
for ((start=1; start<=len; start+=WINDOW_SIZE)); do
    end=$((start + WINDOW_SIZE - 1))
    # If the end position is beyond the sequence length, adjust it
    if ((end > len)); then
        end=$len
    fi

    # Calculate the width (this is just the size of the window, unless it's the last segment)
    width=$((end - start + 1))

    # Output the result to the CSV file
    echo "$id,$seqname,$start,$end,$width" >> $OUTPUT_FILE

    # Increment the ID counter
    ((id++))
done