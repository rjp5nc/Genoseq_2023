#!/bin/bash

# Folder containing .metrics files
METRICS_DIR="/scratch/rjp5nc/UK2022_2024/mapped_bam"
OUTPUT_FILE="collated_metrics.tsv"
FIRST=1

> "$OUTPUT_FILE"  # Clear output file to avoid appending to an old file

for file in "$METRICS_DIR"/*.metrics; do
    # Extract the sample name (filename without extension)
    sample_name=$(basename "$file" .metrics)

    # Process the file, skip metadata lines, and add sample name
    awk -v sample="$sample_name" '
    BEGIN { header_printed = 0 }
    # Skip lines starting with ##
    /^##/ { next }
    # When encountering the first line with data, print header
    NR == 1 && !header_printed {
        print $0 "\tsample"  # Add 'sample' column to the header
        header_printed = 1
    }
    # Print metric data lines and append sample name
    !/^$/ {
        print $0 "\t" sample
    }
    ' "$file" >> "$OUTPUT_FILE"
    
    # After the first file, remove the header to avoid duplication
    if [ $FIRST -eq 1 ]; then
        FIRST=0
    else
        sed -i '1d' "$OUTPUT_FILE"
    fi
done

echo "Collated metrics saved to $OUTPUT_FILE"


#Filter for ___
#Copy
#Add headers back in
