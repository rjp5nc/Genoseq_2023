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


module load bcftools

OUTDIR="/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/gatk_gvcf"
VCF="$OUTDIR/usdobtusa.mito.biallelic.clean.vcf"
MITO="/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/mito_types.csv"

TMPDIR="$OUTDIR/top3_plus_missing_tmp"
mkdir -p "$TMPDIR"

OUTVCF="$OUTDIR/usdobtusa.mito.biallelic.clean.top3_plus_missing.vcf.gz"

# avg DP per sample
bcftools query -f '[%SAMPLE\t%DP\n]' "$VCF" \
| awk '{sum[$1]+=$2; n[$1]++} END{for(s in sum) print s, sum[s]/n[s]}' \
> "$TMPDIR/avgdp.txt"

# mitotype map: sample_clone -> Group
awk -F',' 'NR>1{gsub(/"/,"",$1); gsub(/"/,"",$2); print $1"_clone", $2}' "$MITO" \
| sort -u > "$TMPDIR/mitotype.map"

awk -F',' 'NR>1{print $1, $2}' "$MITO" | sort -u > "$TMPDIR/mitotype.map"




bcftools query -l "$VCF" | sort -u > "$TMPDIR/vcf.samples"

# 2) make list of samples missing mitotype (in VCF but not in mitotype file)
comm -23 "$TMPDIR/vcf.samples" <(cut -d' ' -f1 "$TMPDIR/mitotype.map" | sort -u) \
  > "$TMPDIR/vcf_missing_mitotype.samples"

# 3) top 3 avgDP per mitotype (only among samples that are actually in the VCF)
awk 'NR==FNR{inv[$1]=1; next} inv[$1]{dp[$1]=$2} END{for(s in dp) print s, dp[s]}' \
  "$TMPDIR/vcf.samples" "$TMPDIR/avgdp.txt" > "$TMPDIR/avgdp_in_vcf.txt"

join -1 1 -2 1 \
  <(sort -k1,1 "$TMPDIR/mitotype.map") \
  <(sort -k1,1 "$TMPDIR/avgdp_in_vcf.txt") \
| awk '{print $2, $1, $3}' \
| sort -k1,1 -k3,3nr \
| awk '
  {
    grp=$1
    if (++seen[grp] <= 3) print $2
  }
' > "$TMPDIR/top3_by_mitotype.samples"

# 4) final keep list = (top3 per mitotype) + (missing-mitotype samples)
cat "$TMPDIR/top3_by_mitotype.samples" "$TMPDIR/vcf_missing_mitotype.samples" \
  | sort -u > "$TMPDIR/keep.samples"

echo "Top3-by-mitotype kept: $(wc -l < "$TMPDIR/top3_by_mitotype.samples")"
echo "Missing-mitotype kept: $(wc -l < "$TMPDIR/vcf_missing_mitotype.samples")"
echo "Total kept: $(wc -l < "$TMPDIR/keep.samples")"

bcftools view -S "$TMPDIR/keep.samples" -Oz -o "$OUTDIR/usdobtusa.mito.biallelic.clean.top3_by_mitotype_plus_missing.vcf.gz" "$VCF"
bcftools index -f "$OUTDIR/usdobtusa.mito.biallelic.clean.top3_by_mitotype_plus_missing.vcf.gz"










OUTDIR="/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/gatk_gvcf_top3"
wd="$OUTDIR"   # keep everything in the same output dir

VCF_IN="/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/gatk_gvcf/usdobtusa.mito.biallelic.clean.top3_by_mitotype_plus_missing.vcf.gz"

module load bcftools
module load ruby

cd "$wd"

# ----------------------------
# (0) Index VCF if needed
# ----------------------------
bcftools index -f "$VCF_IN"
tabix -p vcf "$VCF_IN"

# ----------------------------
# (1) Clean again: remove ALT="*"
#     Output must be uncompressed VCF for snapp_prep.rb
# ----------------------------
CLEAN_VCF="$wd/usdobtusa.mito.top3.clean.noSTAR.vcf"

bcftools view \
  -e 'ALT="*"' \
  -Ov \
  -o "$CLEAN_VCF" \
  "$VCF_IN"

# ----------------------------
# (2) Make sample list + 2-col mapping
#     (SnAPP needs "taxon sample_name")
# ----------------------------
SAMPLES_TXT="$wd/snapp_top3.samples.txt"
MAP_2COL="$wd/snapp_top3.samples_2col.txt"

bcftools query -l "$CLEAN_VCF" > "$SAMPLES_TXT"

# Make "sample_clone  sample" mapping (matches what you did before)
awk '{print $1 "_clone", $1}' "$SAMPLES_TXT" > "$MAP_2COL"

# ----------------------------
# (3) Constraints file (single monophyletic NA containing all taxa)
# ----------------------------
CONSTR="$wd/snapp_constraints_top3.txt"

echo -n "monophyletic NA " > "$CONSTR"
sed 's/$/_clone/' "$SAMPLES_TXT" | paste -sd, - >> "$CONSTR"

# ----------------------------
# (4) Run snapp_prep.rb -> SNAPP XML
# ----------------------------
XML_OUT="$wd/snapp.mito.top3.xml"
PREFIX_OUT="$wd/snapp.mito.top3"

ruby /scratch/rjp5nc/snapp5/snapp_prep.rb \
  -v "$CLEAN_VCF" \
  -t "$MAP_2COL" \
  -c "$CONSTR" \
  -x "$XML_OUT" \
  -o "$PREFIX_OUT" \
  -m 10000 \
  -l 10000



# WARNING: The maximum number of SNPs has been set to 10000, which is greater
#     than the number of bi-allelic SNPs with sufficient information (492) for SNAPP.
# WARNING: Excluded 1429 sites with only missing data in one or more species.
# WARNING: Excluded 447 monomorphic sites.

# INFO: Retained 492 bi-allelic sites.