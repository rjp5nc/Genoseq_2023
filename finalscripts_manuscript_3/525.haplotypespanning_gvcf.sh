#ijob -A berglandlab -c10 -p standard --mem=40G


module load bcftools samtools gcc/11.4.0 openmpi/4.1.4 icu R/4.3.1 mafft/7.505

REF="/scratch/rjp5nc/Reference_genomes/mito_reference/usdobtusa_mito.fasta"

# These come from your earlier steps
AVGDP="/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/gatk_gvcf/top3_plus_missing_tmp/avgdp.txt"
MITOMAP="/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/gatk_gvcf/top3_plus_missing_tmp/mitotype.map"

# Directory containing per-sample mito gVCFs
GVCFDIR="/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/gatk_gvcf"

OUTDIR="/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/popart_mafft_gvcf"
mkdir -p "$OUTDIR"/{vcf,consensus,tmp}

# index reference if needed
[ -f "${REF}.fai" ] || samtools faidx "$REF"

# ------------------------------------------------------------
# 1) Build list: highest-depth sample PER mitotype
#    avgdp.txt: "Sample 45.6"   (Sample has NO _clone)
#    mitotype.map: "Sample_clone A"
# ------------------------------------------------------------
TOP_BY_GROUP="$OUTDIR/tmp/top1_per_mitotype.samples"
awk '
  FNR==NR { dp[$1]=$2; next }                          # read avgdp: dp[sample]=depth
  { s=$1; g=$2; sub(/_clone$/,"",s);                   # from mitotype.map: strip _clone
    if (!(s in dp)) next;
    if (!(g in best) || dp[s] > bestdp[g]) { best[g]=s; bestdp[g]=dp[s] }
  }
  END { for (g in best) printf "%s\t%s\t%.6f\n", g, best[g], bestdp[g] }
' "$AVGDP" "$MITOMAP" \
| sort -k1,1 > "$OUTDIR/tmp/top1_per_mitotype.tsv"

cut -f2 "$OUTDIR/tmp/top1_per_mitotype.tsv" > "$TOP_BY_GROUP"

echo "Selected top sample per mitotype:"
cat "$OUTDIR/tmp/top1_per_mitotype.tsv"

# ------------------------------------------------------------
# 2) For each selected sample, build consensus from its gVCF
# ------------------------------------------------------------
ALLFA="$OUTDIR/group_consensus_all.fa"
: > "$ALLFA"

while read -r SAMPLE; do
  # Find gVCF file (adjust patterns if your naming differs)
  G1="${GVCFDIR}/gvcf/${SAMPLE}.g.vcf.gz"
  G2="${GVCFDIR}/gvcf/${SAMPLE}_clone.g.vcf.gz"

  if [ -f "$G1" ]; then
    GVCF="$G1"
    TIP="${SAMPLE}_clone"
  elif [ -f "$G2" ]; then
    GVCF="$G2"
    TIP="${SAMPLE}_clone"
  else
    echo "WARN: missing gVCF for $SAMPLE (tried $G1 / $G2)" >&2
    continue
  fi

  # Make a variant-only VCF (drop gVCF ref blocks TYPE="ref")
  VCF="$OUTDIR/vcf/${SAMPLE}.vcf.gz"
  bcftools view -i 'TYPE!="ref"' -Oz -o "$VCF" "$GVCF"
  bcftools index -f "$VCF"

  # Build consensus (choose haplotype 1 if diploid calls exist; remove -H 1 if you want IUPAC)
  FA="$OUTDIR/consensus/${SAMPLE}.fa"
  bcftools consensus -H 1 -f "$REF" "$VCF" > "$FA"

  # Name sequence by sample (PopART will use this header)
  sed -i "1s/^>.*/>${TIP}/" "$FA"

  cat "$FA" >> "$ALLFA"
  echo "Consensus: ${TIP} from $(basename "$GVCF")"
done < "$TOP_BY_GROUP"

# ------------------------------------------------------------
# 3) Align for PopART
# ------------------------------------------------------------
mafft --auto --thread 8 "$ALLFA" > "$OUTDIR/group_consensus_all.aln.fa"

echo "Wrote:"
echo "  $OUTDIR/tmp/top1_per_mitotype.tsv"
echo "  $ALLFA"
echo "  $OUTDIR/group_consensus_all_gvcf.aln.fa"





R

library(ape)

fa <- read.FASTA("/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/popart_mafft_gvcf/group_consensus_all_gvcf.aln.fa")
write.nexus.data(
  fa,
  file = "/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/popart_mafft_gvcf/group_consensus_all_gvcf.aln.nex",
  format = "dna"
)
