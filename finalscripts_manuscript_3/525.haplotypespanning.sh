#ijob -A berglandlab -c10 -p standard --mem=40G

module load bcftools samtools mafft

REF="/scratch/rjp5nc/Reference_genomes/mito_reference/usdobtusa_mito.fasta"
BAMDIR="/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/final_mitobam_rg3"
TYPES="/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/mito_types.csv"
TOP="/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/top_samples.txt"

OUTDIR="/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/popart_mafft"
mkdir -p "$OUTDIR"/{vcf,consensus}

# index reference if needed
[ -f "${REF}.fai" ] || samtools faidx "$REF"

while read -r SAMPLE; do
  GROUP=$(awk -F',' -v s="$SAMPLE" 'NR>1 && $1==s {print $2; exit}' "$TYPES")
  BAM="${BAMDIR}/${SAMPLE}finalmap_RG.bam"

  if [ -z "$GROUP" ]; then
    echo "WARN: no group for $SAMPLE" >&2
    continue
  fi
  if [ ! -f "$BAM" ]; then
    echo "WARN: missing BAM $BAM" >&2
    continue
  fi

  VCF="$OUTDIR/vcf/${GROUP}.vcf.gz"
  FA="$OUTDIR/consensus/${GROUP}.fa"

  # haploid calls then consensus
  bcftools mpileup -f "$REF" -Ou "$BAM" \
    | bcftools call -c --ploidy 1 -Oz -o "$VCF"
  bcftools index -f "$VCF"
  bcftools consensus -f "$REF" "$VCF" > "$FA"

  # rename header to group label (popART will use these names)
  sed -i "1s/^>.*/>${GROUP}/" "$FA"

  echo "Consensus: $GROUP from $SAMPLE"
done < "$TOP"

cat "$OUTDIR"/consensus/*.fa > "$OUTDIR/group_consensus_all.fa"


mafft --auto --thread 8 "$OUTDIR/group_consensus_all.fa" > "$OUTDIR/group_consensus_all.aln.fa"





R

library(ape)

fa <- read.FASTA("/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/popart_mafft/group_consensus_all.aln.fa")
write.nexus.data(
  fa,
  file = "/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/popart_mafft/group_consensus_all.aln.nex",
  format = "dna"
)
