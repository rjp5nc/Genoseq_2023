cd /scratch/rjp5nc/UK2022_2024/final_mitobam_rg2

module load bcftools

meta="/scratch/rjp5nc/UK2022_2024/touseforDBI_mito_fullref.csv"
reference="/scratch/rjp5nc/Reference_genomes/mito_reference/usdobtusa_mito.fasta"
outdir="/scratch/rjp5nc/UK2022_2024/allsites_mito/vcfs"
merged="/scratch/rjp5nc/UK2022_2024/allsites_mito/usdobtusa_mito_allsites.vcf.gz"

mkdir -p "$outdir"

# Ensure reference is indexed
[ -f "${reference}.fai" ] || samtools faidx "$reference"

# Build BAM list for usdobtusa_mito (skip header row)
awk -F, 'NR>1 {gsub(/[^a-zA-Z0-9_]/,"",$6); if ($6=="usdobtusa_mito") print $1 "finalmap_RG.bam"}' "$meta" \
  > usdobtusa_bams.txt

# Make per-sample haploid VCFs
while read -r bam; do
  [ -s "$bam" ] || { echo "Missing BAM: $bam" >&2; continue; }

  # Prefer sample name from BAM header (SM tag); fallback to filename stem
  sm=$(bcftools query -l "$bam" 2>/dev/null || true)
  if [ -z "$sm" ]; then
    sm=$(basename "$bam" finalmap_RG.bam)
  fi

  # Call haploid; -A keeps all covered sites (remove -A if you only want variant sites)
  bcftools mpileup \
      -f "$reference" \
      -q 20 -Q 20 \
      -a FORMAT/DP,FORMAT/AD \
      -Ou \
      "$bam" \
  | bcftools call \
      -m -A \
      --ploidy 1 \
      -Oz -o "${outdir}/${sm}.allsites.vcf.gz"

  bcftools index -t "${outdir}/${sm}.allsites.vcf.gz"
done < usdobtusa_bams.txt

# Merge per-sample VCFs; maintain haploid GTs
bcftools merge -m none -Oz \
  -o "$merged"  ${outdir}/*.allsites.vcf.gz

bcftools index -t "$merged"

echo "Done. Haploid multi-sample VCF: $merged"
