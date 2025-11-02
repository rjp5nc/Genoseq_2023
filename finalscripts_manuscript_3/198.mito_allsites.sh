cd /scratch/rjp5nc/UK2022_2024/final_mitobam_rg2
module load bcftools samtools

meta="/scratch/rjp5nc/UK2022_2024/touseforDBI_mito_fullref.csv"
reference="/scratch/rjp5nc/Reference_genomes/mito_reference/usdobtusa_mito.fasta"
outdir="/scratch/rjp5nc/UK2022_2024/allsites_mito"
vcf="${outdir}/usdobtusa_mito_allsites.haploid.vcf.gz"

mkdir -p "$outdir"
[ -f "${reference}.fai" ] || samtools faidx "$reference"

# Build BAM list (skip CSV header)
awk -F, 'NR>1 {gsub(/[^a-zA-Z0-9_]/,"",$6); if ($6=="usdobtusa_mito") print $1 "finalmap_RG.bam"}' \
  "$meta" > usdobtusa_bams.txt

# Multi-sample call in one go, haploid, keep all covered sites
bcftools mpileup \
  -f "$reference" \
  -q 20 -Q 20 \
  -a FORMAT/DP,FORMAT/AD \
  -Ou -b usdobtusa_bams.txt \
| bcftools call \
  -m -A \
  --ploidy 1 \
  -Oz -o "$vcf"







bcftools index -t "$vcf"
echo "Wrote haploid all-sites multi-sample VCF -> $vcf"
