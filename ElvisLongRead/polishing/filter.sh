
makeblastdb -in totalHiCwithallbestgapclosed.fa -dbtype nucl -out subject_genome_db
blastn -query dap.contigs.fasta -db subject_genome_db -outfmt 6 -max_target_seqs 1 -out top_hits.txt


awk '!seen[$1]++' top_hits.txt > top_hits_unique.txt


awk '{print $1}' top_hits_unique.txt | seqtk subseq dap.contigs.fasta - > query_hits.fasta


quast.py /scratch/rjp5nc/blastNkeep/query_hits.fasta -o /scratch/rjp5nc/blastNkeep/out
