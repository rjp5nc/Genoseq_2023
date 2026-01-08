#!/usr/bin/env bash

#SBATCH -J windowed_het # A single job name for the array
#SBATCH --ntasks-per-node=2 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 3-0:00:00 ### 15 seconds
#SBATCH --mem 300G
#SBATCH -o /scratch/rjp5nc/erroroutputs/beast.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/erroroutputs/beast.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab
#SBATCH --mail-type=END               # Send email at job completion
#SBATCH --mail-user=rjp5nc@virginia.edu    # Email address for notifications
#SBATCH --array=1-12

module load gatk

REF=/scratch/rjp5nc/Reference_genomes/orig_ref/eu_pulex_totalHiCwithallbestgapclosed.clean.fa

cd /scratch/rjp5nc/UK2022_2024/daphnia_phylo/eupulex_indv/gvcf/

# ls /scratch/rjp5nc/UK2022_2024/daphnia_phylo/gvcf/eupulex_chr \
#   | grep '^Scaffold_' \
#   | sort -V > scaffolds.txt


SCAF=$(sed -n "${SLURM_ARRAY_TASK_ID}p" scaffolds.txt)

IN_DIR=/scratch/rjp5nc/UK2022_2024/daphnia_phylo/gvcf/eupulex_chr/${SCAF}
OUT_DIR=/scratch/rjp5nc/UK2022_2024/daphnia_phylo/eupulex_indv/gvcf/genomicsdb/${SCAF}
mkdir -p "$OUT_DIR"

# list inputs for this scaffold
find "$IN_DIR" -name "*.g.vcf.gz" > ${SCAF}.gvcfs.list

# (Important) ensure indexes exist, otherwise everything crawls
# ls *.tbi should exist for each g.vcf.gz; if not, create once:
# while read f; do tabix -p vcf "$f"; done < ${SCAF}.gvcfs.list

gatk --java-options "-Xmx56g -Djava.io.tmpdir=$OUT_DIR/tmp" GenomicsDBImport \
  -R "$REF" \
  --genomicsdb-workspace-path "$OUT_DIR/db" \
  --sample-name-map <(awk -v d="$IN_DIR" '
      BEGIN{OFS="\t"}
      { 
        g=$0
        # sample name from filename; adjust if your naming differs
        n=g; sub(/^.*\//,"",n); sub(/\.g\.vcf\.gz$/,"",n)
        print n, g
      }' ${SCAF}.gvcfs.list) \
  -L "$SCAF" \
  --reader-threads 8 \
  --batch-size 50

gatk --java-options "-Xmx56g -Djava.io.tmpdir=$OUT_DIR/tmp" GenotypeGVCFs \
  -R "$REF" \
  -V "gendb://$OUT_DIR/db" \
  -L "$SCAF" \
  -O ${SCAF}.vcf.gz
