#!/usr/bin/env bash
#
#SBATCH -J genotypegvcf # A single job name for the array
#SBATCH --ntasks-per-node=10 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 7:00:00 # 8 hours
#SBATCH --mem 80G
#SBATCH -o /scratch/rjp5nc/err/genotypegvcf.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/err/genotypegvcf.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab
#SBATCH --array=1-112%50
#SBATCH --mail-type=END               # Send email at job completion
#SBATCH --mail-user=rjp5nc@virginia.edu    # Email address for notifications

module load gatk/4.6.0.0
module load samtools/1.17

REF="/scratch/rjp5nc/Reference_genomes/mito_reference/eudobtusa_mito_reverse.fasta"
GVCFDIR="/scratch/rjp5nc/UK2022_2024/redone_mito/euobtusa/gvcf_obtusa"
OUTDIR="/scratch/rjp5nc/UK2022_2024/redone_mito/euobtusa/cohort_gendb"
mkdir -p "$OUTDIR"

SAMPLEMAP="${OUTDIR}/sample_map.tsv"

ls -1 "${GVCFDIR}"/*.g.vcf.gz | sort \
| awk -F'/' '{
    f=$NF
    sub(/\.g\.vcf\.gz$/,"",f)
    print f "\t" $0
  }' > "$SAMPLEMAP"

echo "Wrote: $SAMPLEMAP"
wc -l "$SAMPLEMAP"
head "$SAMPLEMAP"


awk -F'\t' '{print $2}' "$SAMPLEMAP" | while read -r f; do
  [[ -s "${f}.tbi" ]] || echo "MISSING_TBI $f"
done | head


module load samtools/1.17
samtools faidx "$REF"

ALLBED="${OUTDIR}/all_sites.bed"
cut -f1,2 "${REF}.fai" | awk '{print $1"\t0\t"$2}' > "$ALLBED"
echo "Wrote: $ALLBED"
head "$ALLBED"


module load gatk/4.6.0.0
module load samtools/1.17

THREADS="${SLURM_CPUS_PER_TASK:-1}"

REF="/scratch/rjp5nc/Reference_genomes/mito_reference/eudobtusa_mito_reverse.fasta"
OUTDIR="/scratch/rjp5nc/UK2022_2024/redone_mito/euobtusa/cohort_gendb"
SAMPLEMAP="${OUTDIR}/sample_map.tsv"
ALLBED="${OUTDIR}/all_sites.bed"
DB="${OUTDIR}/gendb_mito_obtusa"

mkdir -p "$OUTDIR"

# ensure dict + fai
[[ -s "${REF}.fai" ]] || samtools faidx "$REF"
DICT="${REF%.*}.dict"
[[ -s "$DICT" ]] || gatk CreateSequenceDictionary -R "$REF"

[[ -s "$SAMPLEMAP" ]] || { echo "ERROR: missing $SAMPLEMAP"; exit 1; }
[[ -s "$ALLBED" ]] || { echo "ERROR: missing $ALLBED"; exit 1; }

echo "[$(date)] Importing to GenomicsDB: $DB"
echo "Samples: $(wc -l < "$SAMPLEMAP")"
echo "Intervals: $ALLBED"

gatk --java-options "-Xmx56g" GenomicsDBImport \
  --genomicsdb-workspace-path "$DB" \
  --sample-name-map "$SAMPLEMAP" \
  -L "$ALLBED" \
  --reader-threads "$THREADS" \
  --batch-size 50

echo "[$(date)] Done GenomicsDBImport"









REF="/scratch/rjp5nc/Reference_genomes/mito_reference/eudobtusa_mito_reverse.fasta"
OUTDIR="/scratch/rjp5nc/UK2022_2024/redone_mito/euobtusa/cohort_gendb"
ALLBED="${OUTDIR}/all_sites.bed"
DB="${OUTDIR}/gendb_mito_obtusa"

VCF="${OUTDIR}/obtusa.mito.ALLSITES.vcf.gz"

# ensure dict + fai
[[ -s "${REF}.fai" ]] || samtools faidx "$REF"
DICT="${REF%.*}.dict"
[[ -s "$DICT" ]] || gatk CreateSequenceDictionary -R "$REF"

[[ -d "$DB" ]] || { echo "ERROR: missing GenomicsDB workspace: $DB"; exit 1; }
[[ -s "$ALLBED" ]] || { echo "ERROR: missing $ALLBED"; exit 1; }

echo "[$(date)] Genotyping ALL SITES -> $VCF"

gatk --java-options "-Xmx40g" GenotypeGVCFs \
  -R "$REF" \
  -V "gendb://$DB" \
  -L "$ALLBED" \
  --include-non-variant-sites \
  -O "$VCF"

echo "[$(date)] Done GenotypeGVCFs"
ls -lh "$VCF"*






module load bcftools || true

REFLEN=$(awk '{s+=$2} END{print s}' "${REF}.fai")
echo "Reference length sum: $REFLEN"

# count VCF records (non-header)
zcat /scratch/rjp5nc/UK2022_2024/redone_mito/euobtusa/cohort_gendb/obtusa.mito.ALLSITES.vcf.gz \
| grep -v '^#' | wc -l
