#!/usr/bin/env bash
#
#SBATCH -J run_blast # A single job name for the array
#SBATCH --ntasks-per-node=10 # multi core
#SBATCH -N 1 # on one node
#SBATCH -t 0-03:00 # 3 hours
#SBATCH --mem 40G
#SBATCH -o /scratch/rjp5nc/Canu_error/blast%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/Canu_error/blast%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# Load required modules (if needed)
module load blast

# Define directories
INPUT_DIR="/scratch/rjp5nc/HMW/HMWDNAElvis3/blast/trimmedfasta"
OUTPUT_DIR="/scratch/rjp5nc/HMW/HMWDNAElvis3/blast_results2"
mkdir -p "$OUTPUT_DIR"

# Get list of FASTA files
FASTA_FILES=($INPUT_DIR/*.fasta)

# Get the specific file for this array job
FASTA_FILE=${FASTA_FILES[$SLURM_ARRAY_TASK_ID-1]}
BASENAME=$(basename "$FASTA_FILE" .fasta)

# Define output files
BLAST_OUTPUT="$OUTPUT_DIR/${BASENAME}_blast.txt"
ACCESSION_FILE="$OUTPUT_DIR/${BASENAME}_accessions.txt"
FINAL_OUTPUT="$OUTPUT_DIR/${BASENAME}_final_results.txt"

echo "Processing: $FASTA_FILE"
echo "Saving BLAST results to: $BLAST_OUTPUT"

# Run BLASTn remotely
blastn -query "$FASTA_FILE" -db "nt" -remote -evalue 1e-5 -max_target_seqs 1 \
    -outfmt "6 qacc sacc pident qstart qend sstart send" \
    -out "$BLAST_OUTPUT"

# Check if BLAST produced any output
if [[ ! -s "$BLAST_OUTPUT" ]]; then
    echo "BLAST failed or returned no results for $FASTA_FILE"
    exit 1
fi

# Extract query accession (qacc) and subject accession (sacc)
awk '{print $1, $2}' "$BLAST_OUTPUT" | sort -u > "$ACCESSION_FILE"

echo "Retrieving species names for accession numbers..."

# Process each accession number
> "$FINAL_OUTPUT"
while read -r qacc accession; do
    echo "Fetching taxonomy for: $accession"
    
    species=$(esearch -db nucleotide -query "$accession" | efetch -format docsum | xtract -pattern DocumentSummary -element Title)
    
    if [[ -z "$species" ]]; then
        species="Unknown"
    fi
    
    echo -e "$qacc\t$accession\t$species" >> "$FINAL_OUTPUT"
done < "$ACCESSION_FILE"

echo "✅ Completed processing $FASTA_FILE. Results saved to: $FINAL_OUTPUT"
