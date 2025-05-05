module load samtools

for bam in /scratch/rjp5nc/UK2022_2024/allshortreads/alignedbamsnew/*.bam; do
    samtools index -@ 10 "$bam" "${bam}.bai"
done


for bam in /scratch/rjp5nc/UK2022_2024/allshortreads/alignedbamsnew/*.bam; do
    samtools quickcheck "$bam" || echo "$bam has issues"
done