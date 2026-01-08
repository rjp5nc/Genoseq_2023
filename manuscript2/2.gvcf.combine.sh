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

cd /scratch/rjp5nc/UK2022_2024/daphnia_phylo/eupulex_indv/gvcf/


module load gatk

REF=/scratch/rjp5nc/Reference_genomes/orig_ref/eu_pulex_totalHiCwithallbestgapclosed.clean.fa
SCAF=$(sed -n "${SLURM_ARRAY_TASK_ID}p" scaffolds.txt)

IN_DIR=/scratch/rjp5nc/UK2022_2024/daphnia_phylo/gvcf/eupulex_chr/${SCAF}
OUT_DIR=/scratch/rjp5nc/UK2022_2024/daphnia_phylo/eupulex_indv/gvcf/combined_by_scaffold
TMP_DIR=/scratch/rjp5nc/UK2022_2024/daphnia_phylo/eupulex_indv/gvcf/tmp/${SCAF}

mkdir -p "$OUT_DIR" "$TMP_DIR"

find "$IN_DIR" -name "*.g.vcf.gz" | sort > "$TMP_DIR/${SCAF}.gvcfs.list"

# optional but recommended: check indexes
while read f; do
  [[ -f "${f}.tbi" ]] || { echo "Missing index: ${f}.tbi" >&2; exit 2; }
done < "$TMP_DIR/${SCAF}.gvcfs.list"

# CombineGVCFs with xargs to avoid huge command line expansions
xargs -a "$TMP_DIR/${SCAF}.gvcfs.list" -I{} echo -V {} > "$TMP_DIR/${SCAF}.V.args"

gatk --java-options "-Xmx56g -Djava.io.tmpdir=$TMP_DIR" CombineGVCFs \
  -R "$REF" \
  $(cat "$TMP_DIR/${SCAF}.V.args") \
  -L "$SCAF" \
  -O "$OUT_DIR/${SCAF}.merged.g.vcf.gz"
