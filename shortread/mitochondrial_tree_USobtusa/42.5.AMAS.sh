
ijob -A berglandlab -c10 -p standard --mem=40G

conda create -n amas-env python=3.11 -y
conda activate amas-env

/scratch/rjp5nc/AMAS

python /scratch/rjp5nc/AMAS/amas/AMAS.py summary -f fasta -d dna \
-i /scratch/rjp5nc/UK2022_2024/consensusmitoaligned2/all_aligned_unique.fasta \
-o /scratch/rjp5nc/UK2022_2024/mito_summary.txt