
#us obtusa, eu obtusa, us ambigua, eu pulex, eu pulex2

us_pulex_ref_kap4.fa.out

tail -n +4 /scratch/rjp5nc/Reference_genomes/masked/Daphnia_ambigua_Q001_genome.fa.out \
| awk 'BEGIN{OFS= "\t"} {print $5, $6, $7}' > /scratch/rjp5nc/UK2022_2024/daphnia_phylo/mask_beds/us_ambigua_mask.bed



output=trimmed10bp_masked_usambigua.vcf

bedtools subtract \
-header \
-a /scratch/rjp5nc/UK2022_2024/daphnia_phylo/trimmed10bp_vcf/trimmed10bp_usambigua_vcf.vcf.gz \
-b /scratch/rjp5nc/UK2022_2024/daphnia_phylo/mask_beds/us_ambigua_mask.bed \
> /scratch/rjp5nc/UK2022_2024/daphnia_phylo/trimmed_10bp_repeatmasked_vcf/$output


bgzip /scratch/rjp5nc/UK2022_2024/daphnia_phylo/trimmed_10bp_repeatmasked_vcf/$output

module load bcftools

bcftools index /scratch/rjp5nc/UK2022_2024/daphnia_phylo/trimmed_10bp_repeatmasked_vcf/$output.gz







output=trimmed10bp_masked_usobtusa.vcf

samples=($(bcftools query -l /scratch/rjp5nc/UK2022_2024/daphnia_phylo/trimmed_10bp_repeatmasked_vcf/$output.gz))

# Build the environment string for sample name mappings
env_vars=$(for i in "${!samples[@]}"; do echo "SAMPLE_$((i+1))=${samples[$i]}"; done)

# Run bcftools and awk with env
bcftools query -f '%CHROM\t%POS[\t%DP]\n' /scratch/rjp5nc/UK2022_2024/daphnia_phylo/trimmed_10bp_repeatmasked_vcf/$output.gz | \
  env $env_vars awk -v OFS="\t" -v n="${#samples[@]}" '
  {
    for(i=3; i<=NF; i++) {
      sum[i-2] += $i; count[i-2]++;
    }
  }
  END {
    for(i=1; i<=n; i++) {
      sample = ENVIRON["SAMPLE_" i];
      avg = (count[i] > 0) ? sum[i] / count[i] : 0;
      printf "%s\t%.2f\n", sample, avg;
    }
  }' > /scratch/rjp5nc/UK2022_2024/daphnia_phylo/trimmed_10bp_repeatmasked_vcf/${output}_depth.txt



#per basepair

bcftools query -f '%CHROM\t%POS[\t%DP]\n' /scratch/rjp5nc/UK2022_2024/daphnia_phylo/trimmed_10bp_repeatmasked_vcf/$output.gz | \
awk '
{
  sum = 0; count = 0;
  for(i=3; i<=NF; i++) {
    if ($i != ".") { sum += $i; count++ }
  }
  avg = (count > 0) ? sum / count : 0;
  printf "%s\t%s\t%.2f\n", $1, $2, avg;
}' > /scratch/rjp5nc/UK2022_2024/daphnia_phylo/trimmed_10bp_repeatmasked_vcf/${output}_depth_per_site.txt