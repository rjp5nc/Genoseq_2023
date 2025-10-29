cd /scratch/rjp5nc/UK2022_2024/final_mitobam_rg2

module load bcftools

meta="/scratch/rjp5nc/UK2022_2024/touseforDBI_mito_fullref.csv"

# Extract BAMs corresponding to usdobtusa_mito
awk -F, '{gsub(/[^a-zA-Z0-9_]/,"",$6); if($6=="usdobtusa_mito") print $1"finalmap_RG.bam"}' /scratch/rjp5nc/UK2022_2024/touseforDBI_mito_fullref.csv > usdobtusa_bams.txt

# Create all-sites VCFs for each BAM
reference="/scratch/rjp5nc/Reference_genomes/mito_reference/usdobtusa_mito.fasta"

while read bam; do
    sample=$(basename $bam finalmap_RG.bam)
    bcftools mpileup \
        -f $reference \
        -a FORMAT/DP,AD \
        -Ou \
        $bam | \
    bcftools call \
        -mv \
        --ploidy 1 \
        -Oz \
        -o /scratch/rjp5nc/UK2022_2024/allsites_mito/vcfs/${sample}.allsites.vcf.gz

    # Index the output
    bcftools index -t /scratch/rjp5nc/UK2022_2024/allsites_mito/vcfs/${sample}.allsites.vcf.gz
done < usdobtusa_bams.txt

cd /scratch/rjp5nc/UK2022_2024/allsites_mito/vcfs/

bcftools merge -m none -Oz -o /scratch/rjp5nc/UK2022_2024/allsites_mito/usdobtusa_mito_allsites.vcf.gz  /scratch/rjp5nc/UK2022_2024/allsites_mito/vcfs/*.allsites.vcf.gz


# Index the merged VCF
bcftools index -t /scratch/rjp5nc/UK2022_2024/allsites_mito/usdobtusa_mito_allsites.vcf.gz