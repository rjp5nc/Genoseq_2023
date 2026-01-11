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


VCF="/scratch/rjp5nc/UK2022_2024/redone_mito/euobtusa/cohort_gendb/obtusa.mito.ALLSITES.vcf.gz"
KEEP="/scratch/rjp5nc/UK2022_2024/redone_mito/euobtusa/samples_lt30pct_missing.txt"
DROP="/scratch/rjp5nc/UK2022_2024/redone_mito/euobtusa/samples_ge30pct_missing.txt"

tmpdir=$(mktemp -d)
trap 'rm -rf "$tmpdir"' EXIT

cd /scratch/rjp5nc/UK2022_2024/redone_mito/euobtusa/cohort_gendb

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
