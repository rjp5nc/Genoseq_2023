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


WINDOW_SIZE=1000000

echo "id,seqnames,start,end,width" > "$OUTPUT_FILE"
id=1
seqname=""
seq=""

while read -r line; do
    if [[ $line == ">"* ]]; then
        # Process previous scaffold if present
        if [[ -n $seq ]]; then
            len=${#seq}
            for ((start=1; start<=len; start+=WINDOW_SIZE)); do
                end=$((start + WINDOW_SIZE - 1))
                (( end > len )) && end=$len
                width=$((end - start + 1))
                echo "$id,$seqname,$start,$end,$width" >> "$OUTPUT_FILE"
                ((id++))
            done
        fi
        # Start new scaffold
        seqname=${line#>}
        seq=""
    else
        seq+=$line
    fi
done < "$FASTA_FILE"

# Final scaffold processing
if [[ -n $seq ]]; then
    len=${#seq}
    for ((start=1; start<=len; start+=WINDOW_SIZE)); do
        end=$((start + WINDOW_SIZE - 1))
        (( end > len )) && end=$len
        width=$((end - start + 1))
        echo "$id,$seqname,$start,$end,$width" >> "$OUTPUT_FILE"
        ((id++))
    done
fi
