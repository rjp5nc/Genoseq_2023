
cd /scratch/rjp5nc/UK2022_2024/daphnia_phylo/gvcf/euobtusa_chr/h2tg000002l_2

#!/bin/bash

# Fix all .g.vcf.gz files under current directory and subdirectories
#sed to remove the path  : /scratch/rjp5nc/UK2022_2024/mapped_bam/

find . -name "*.g.vcf.gz" | while read gzfile; do
    echo "Processing $gzfile"
    
    # Uncompress to temporary file
    gunzip -c file.g.vcf.gz | sed 's/\r$//' > file.g.vcf
    bgzip -f file.g.vcf
       
    # Remove paths and suffixes from within the file content
    sed -i 's|/scratch/rjp5nc/UK2022_2024/mapped_bam/||g; s|_finalmap.bam||g' "${gzfile%.gz}"
    
    gatk IndexFeatureFile -I file.g.vcf.gz
    # Recompress back
    bgzip -f "${gzfile%.gz}"
done