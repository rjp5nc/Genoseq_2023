#!/usr/bin/env bash

#SBATCH -J windowed_het # A single job name for the array
#SBATCH --ntasks-per-node=10 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 3-0:00:00 ### 15 seconds
#SBATCH --mem 100G
#SBATCH -o /scratch/rjp5nc/erroroutputs/beast.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/erroroutputs/beast.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab
#SBATCH --mail-type=END               # Send email at job completion
#SBATCH --mail-user=rjp5nc@virginia.edu    # Email address for notifications
#SBATCH --array=1-12

# ls /scratch/rjp5nc/UK2022_2024/daphnia_phylo/gvcf/eupulex_chr \
#   | grep '^Scaffold_' \
#   | sort -V > scaffolds.txt


module load gatk
cd /scratch/rjp5nc/UK2022_2024/daphnia_phylo/eupulex_indv/gvcf/

REF=/scratch/rjp5nc/Reference_genomes/orig_ref/eu_pulex_totalHiCwithallbestgapclosed.clean.fa
SCAF=$(sed -n "${SLURM_ARRAY_TASK_ID}p" scaffolds.txt)

IN_DIR=/scratch/rjp5nc/UK2022_2024/daphnia_phylo/gvcf/eupulex_chr/${SCAF}
OUT_DIR=/scratch/rjp5nc/UK2022_2024/daphnia_phylo/eupulex_indv/gvcf/genomicsdb/${SCAF}

mkdir -p "$OUT_DIR" "$OUT_DIR/tmp"

# list gvcfs for this scaffold
find "$IN_DIR" -name "*.g.vcf.gz" | sort > "$OUT_DIR/${SCAF}.gvcfs.list"

# build sample map (tab-separated)
awk -F/ 'BEGIN{OFS="\t"}
{
  fn=$NF
  sub(/\.g\.vcf\.gz$/,"",fn)
  print fn, $0
}' "$OUT_DIR/${SCAF}.gvcfs.list" > "$OUT_DIR/${SCAF}.sample_map.tsv"

# optional: fail early if sample names collide
dups=$(cut -f1 "$OUT_DIR/${SCAF}.sample_map.tsv" | sort | uniq -d | head -n 1)
if [[ -n "$dups" ]]; then
  echo "ERROR: duplicate sample name(s) found for $SCAF (example: $dups)" >&2
  exit 2
fi

# IMPORTANT: remove any previous partial workspace
rm -rf "$OUT_DIR/db"

# run import once
if gatk --java-options "-Xmx56g -Djava.io.tmpdir=$OUT_DIR/tmp" GenomicsDBImport \
    -R "$REF" \
    --genomicsdb-workspace-path "$OUT_DIR/db" \
    --sample-name-map "$OUT_DIR/${SCAF}.sample_map.tsv" \
    -L "$SCAF" \
    --reader-threads 8 \
    --batch-size 50
then
  gatk --java-options "-Xmx56g -Djava.io.tmpdir=$OUT_DIR/tmp" GenotypeGVCFs \
    -R "$REF" \
    -V "gendb://$OUT_DIR/db" \
    -L "$SCAF" \
    -O "$OUT_DIR/${SCAF}.vcf.gz"
else
  echo "GenomicsDBImport failed for $SCAF" >&2
  exit 1
fi
