#!/bin/bash

conda activate blast
module load blast 

#create input fastas
awk '/^>/ {filename = "output_fasta/" substr($1, 2) ".fasta"; print $0 > filename; next} {print $0 > filename}' dap.contigs.fasta




# Define directories
FASTA_DIR="/scratch/rjp5nc/HMW/HMWDNAElvis3/blast/output_fasta"
OUTPUT_DIR="/scratch/rjp5nc/HMW/HMWDNAElvis3/blast/trimmedfasta"
mkdir -p "$OUTPUT_DIR"

# Loop through each .fasta file in the directory
for FASTA_FILE in "$FASTA_DIR"/*.fasta; do
    # Get the length of the sequence in the file
    SEQUENCE_LENGTH=$(grep -v "^>" "$FASTA_FILE" | tr -d '\n' | wc -c)
    
    # Ensure the sequence is long enough (at least 100 bp)
    if [ "$SEQUENCE_LENGTH" -ge 100 ]; then
        # Calculate the middle of the sequence
        MID_POINT=$(( SEQUENCE_LENGTH / 2 ))

        # Choose a random position for the start of the fragment
        START_POS=$(( MID_POINT - 5000 ))

        # Adjust start position if sequence is too short
        if [ "$START_POS" -lt 0 ]; then
            START_POS=0
        fi
        if [ "$START_POS" -gt $(( SEQUENCE_LENGTH - 10000 )) ]; then
            START_POS=$(( SEQUENCE_LENGTH - 10000 ))
        fi

        # Extract the sequence name
        SEQUENCE_NAME=$(head -n 1 "$FASTA_FILE" | sed 's/^>//')
        OUTPUT_FILE="$OUTPUT_DIR/${SEQUENCE_NAME}_random_10k.fasta"
        
        # Extract the sequence fragment and save it
        awk -v start="$START_POS" -v len=10000 -v seq_len="$SEQUENCE_LENGTH" '
        BEGIN { FS = "\n"; RS = ">" }
        NR > 1 {
            seq = $2;
            print ">"$1;
            print substr(seq, start+1, (seq_len < 10000 ? seq_len : len));
        }
        ' "$FASTA_FILE" > "$OUTPUT_FILE"
    fi
done

# Run BLAST on extracted sequences






