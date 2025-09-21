

#do analyzekraken.r first

sed 's/^>//' /scratch/rjp5nc/krakenDB/US_obtusa/US_obtusa_dap_contigs.txt > /scratch/rjp5nc/krakenDB/US_obtusa/US_obtusa_dap_contigs_cleaned.txt


FASTA_FILE="/scratch/rjp5nc/Reference_genomes/orig_ref/JAACYE01.1.fasta"
OUTPUT_FILE="/scratch/rjp5nc/Reference_genomes/post_kraken/US_obtusa_onlydaps.fa"
CONTIG_LIST="/scratch/rjp5nc/krakenDB/US_obtusa/US_obtusa_dap_contigs_cleaned.txt"

# Extract contigs based on the list of names in hap2_daps_cleaned.txt
awk 'BEGIN {while(getline < "'$CONTIG_LIST'") list[$1] = 1} 
     /^>/ {header = substr($0, 2)} 
     header in list {print $0; print_seq=1; next} 
     print_seq == 1 && !/^>/ {print $0} 
     /^>/ {print_seq=0}' $FASTA_FILE > $OUTPUT_FILE


