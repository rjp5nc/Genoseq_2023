#!/usr/bin/env bash
#SBATCH -J DownloadMap    # Job name
#SBATCH --ntasks=1        # Single task per job
#SBATCH --cpus-per-task=10 # Number of CPU cores per task
#SBATCH -N 1              # Run on one node
#SBATCH -t 1:20:00        # 10 hours runtime
#SBATCH --mem=20G        # Memory per node
#SBATCH -o /scratch/rjp5nc/outputerrors/down.%A_%a.out  # Standard output
#SBATCH -e /scratch/rjp5nc/outputerrors/down.%A_%a.err  # Standard error
#SBATCH -p standard       # Partition
#SBATCH --account=berglandlab


VCF="/scratch/rjp5nc/UK2022_2024/redone_mito/eupulex/cohort_gendb/pulex.mito.ALLSITES.vcf.gz"
KEEP="/scratch/rjp5nc/UK2022_2024/redone_mito/eupulex/samples_lt30pct_missing.txt"
DROP="/scratch/rjp5nc/UK2022_2024/redone_mito/eupulex/samples_ge30pct_missing.txt"

tmpdir=$(mktemp -d)
trap 'rm -rf "$tmpdir"' EXIT

cd /scratch/rjp5nc/UK2022_2024/redone_mito/eupulex/cohort_gendb

# list samples
bcftools query -l "$VCF" > "$tmpdir/samples.txt"

# compute missing fraction for each sample (./. , .|. , . treated as missing)
: > "$KEEP"
: > "$DROP"

while read -r S; do
  frac=$(bcftools view -s "$S" "$VCF" -Ou \
    | bcftools query -f '[%GT\n]' \
    | awk '
        BEGIN{n=0;m=0}
        { n++; g=$1; if(g=="./." || g==".|." || g==".") m++ }
        END{ if(n==0) print "NA"; else printf "%.6f", m/n }')

  # keep if frac is numeric and <0.30
  awk -v s="$S" -v f="$frac" 'BEGIN{
      if(f=="NA"){exit}
      if(f+0 < 0.30) print s "\t" f;
      else           print s "\t" f > "/dev/stderr"
    }' >> "$KEEP" 2>> "$DROP"
done < "$tmpdir/samples.txt"

# Keep/drop are currently "sample<TAB>frac" for transparency
echo "Keep (<30% missing): $(wc -l < "$KEEP")"
echo "Drop (>=30% missing or NA): $(wc -l < "$DROP")"

# If you want KEEP to be sample names only:
cut -f1 "$KEEP" > "${KEEP%.txt}.names.txt"
echo "Wrote: ${KEEP%.txt}.names.txt"



bcftools view -H /scratch/rjp5nc/UK2022_2024/redone_mito/eupulex/cohort_gendb/pulex.mito.ALLSITES.vcf.gz | head


bcftools query -f '[%GT\t]\n' pulex.mito.ALLSITES.vcf.gz | \
awk '{
  for(i=1;i<=NF;i++){
    if($i=="./.") miss[i]++
    tot[i]++
  }
}
END{
  for(i=1;i<=NF;i++) printf "sample_col_%d\t%d\t%d\t%.4f\n", i, miss[i], tot[i], miss[i]/tot[i]
}' > nocall_by_column.txt


VCF=/scratch/rjp5nc/UK2022_2024/redone_mito/eupulex/cohort_gendb/pulex.mito.ALLSITES.vcf.gz

# sample list in order (1-based like your sample_col_1,2,...)
bcftools query -l "$VCF" | nl -w1 -s$'\t' > sample_order.txt

# join with your nocall_by_column.txt
awk 'BEGIN{OFS="\t"} {sub("sample_col_","",$1); print $1,$2,$3,$4}' nocall_by_column.txt > nocall_tmp.txt

join -t $'\t' -1 1 -2 1 sample_order.txt nocall_tmp.txt | \
awk 'BEGIN{OFS="\t"} {print $2,$1,$3,$4,$5}' \
> nocall_by_sample.txt

# Columns: sample_name  col_index  n_nocall  n_total  frac_nocall
sort -k5,5nr nocall_by_sample.txt | head -50




awk '$5 > 0.3 {print $1}' nocall_by_sample.txt | sort -u > remove_samples_gt0.3.txt
wc -l remove_samples_gt0.3.txt
head remove_samples_gt0.3.txt



VCF=/scratch/rjp5nc/UK2022_2024/redone_mito/eupulex/cohort_gendb/pulex.mito.ALLSITES.vcf.gz
OUT=/scratch/rjp5nc/UK2022_2024/redone_mito/eupulex/cohort_gendb/pulex.mito.ALLSITES.rmMissGT0.3.vcf.gz

bcftools view -S ^remove_samples_gt0.3.txt -Oz -o "$OUT" "$VCF"
bcftools index -t "$OUT"
