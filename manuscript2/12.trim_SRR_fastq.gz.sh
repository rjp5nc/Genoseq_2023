#!/usr/bin/env bash
#
#SBATCH -J trimSRA # A single job name for the array
#SBATCH --cpus-per-task=10 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-10:00:00 # 10 hours
#SBATCH --mem 100G
#SBATCH -o /scratch/rjp5nc/erroroutputs/down.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/erroroutputs/down.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab
#SBATCH --array=1-643%60   # Adjust the range based on the number of folders


#mv /scratch/rjp5nc/UK2022_2024/allshortreads/01.RawData/SR* /scratch/rjp5nc/UK2022_2024/allshortreads/01.RawData/SRR/



module load cutadapt gcc/11.4.0 bwa/0.7.17 samtools/1.17 picard/2.27.5

# -----------------------
# CONFIG
# -----------------------
#IDFILE="/scratch/rjp5nc/rawdata/sra_ids.txt"
IDFILE="/scratch/rjp5nc/rawdata/sra_metadata_out/sra_ids_dobtusa.txt"
RAWBASE="/scratch/rjp5nc/rawdata/SRRsamps"
REF="/scratch/rjp5nc/Reference_genomes/mito_reference/eudobtusa_mito_reverse.fasta"


# cd /scratch/rjp5nc/rawdata/sra_metadata_out

# awk -F'\t' '
#   NR>1 &&
#   $1 ~ /^SRR/ &&
#   $6 == "Daphnia obtusa"
#   { print $1 }
# ' sra_merged.tsv | sort -u > ../sra_ids_dobtusa.txt



# cd /scratch/rjp5nc/rawdata/sra_metadata_out

# awk -F'\t' '
#   NR>1 &&
#   $1 ~ /^SRR/ &&
#   $6 == "Daphnia pulex"
#   { print $1 }
# ' sra_merged.tsv | sort -u > ../sra_ids_dpulex.txt



OUTBASE="/scratch/rjp5nc/UK2022_2024/redone_mito/euobtusa/SRR"
mkdir -p "${OUTBASE}/trimmed_fastq" \
         "${OUTBASE}/sortedbams" \
         "${OUTBASE}/sortedbamsdedup" \
         "${OUTBASE}/sortedbamsreport"

# -----------------------
# GET SRR FOR THIS TASK
# -----------------------
samp=$(awk -v n="${SLURM_ARRAY_TASK_ID}" 'NR==n{gsub(/\r/,""); print $1; exit}' "${IDFILE}")

if [[ -z "${samp}" ]]; then
  echo "Error: empty SRR at task ${SLURM_ARRAY_TASK_ID}"
  exit 1
fi

dir="${RAWBASE}/${samp}"
R1="${dir}/${samp}_1.fastq.gz"
R2="${dir}/${samp}_2.fastq.gz"

echo "[$(date)] Task ${SLURM_ARRAY_TASK_ID}: ${samp}"
echo "Dir: ${dir}"
echo "R1:  ${R1}"
echo "R2:  ${R2}"

if [[ ! -s "${R1}" || ! -s "${R2}" ]]; then
  echo "Error: missing input FASTQs for ${samp}"
  ls -lh "${dir}" || true
  exit 1
fi

# -----------------------
# TRIM (cutadapt)
# -----------------------
trim1="${OUTBASE}/trimmed_fastq/${samp}.trimmed1.fq.gz"
trim2="${OUTBASE}/trimmed_fastq/${samp}.trimmed2.fq.gz"

cutadapt \
  -q 18 \
  --minimum-length 75 \
  -O 15 \
  -n 3 \
  --cores="${SLURM_CPUS_PER_TASK}" \
  -o "${trim1}" \
  -p "${trim2}" \
  "${R1}" "${R2}"

# -----------------------
# MAP + SORT
# -----------------------
sortedbam="${OUTBASE}/sortedbams/${samp}.sorted.bam"

bwa mem \
  -t "${SLURM_CPUS_PER_TASK}" \
  -R "@RG\tID:${samp}\tSM:${samp}\tPL:illumina\tLB:lib1" \
  "${REF}" \
  "${trim1}" "${trim2}" | \
samtools view -@ "${SLURM_CPUS_PER_TASK}" -Sbh -q 20 -F 0x100 - | \
samtools sort -@ "${SLURM_CPUS_PER_TASK}" -o "${sortedbam}" -

samtools index "${sortedbam}"

# -----------------------
# MARK DUPLICATES (remove PCR dups)
# -----------------------
dedupbam="/scratch/rjp5nc/UK2022_2024/redone_mito/euobtusa/sortedbamsdedup_obtusa/${samp}.sort.dedup.bam"
metrics="${OUTBASE}/sortedbamsreport/${samp}.mark_duplicates_report.txt"

java -jar "$EBROOTPICARD/picard.jar" MarkDuplicates \
  -REMOVE_DUPLICATES true \
  -I "${sortedbam}" \
  -O "${dedupbam}" \
  -M "${metrics}" \
  -VALIDATION_STRINGENCY SILENT \
  -CREATE_INDEX true

echo "[$(date)] Finished ${samp}"


# samtools stats /scratch/rjp5nc/UK2022_2024/redone_mito/euobtusa/sortedbamsdedup/${samp}.sort.dedup.bam > /scratch/rjp5nc/"${samp}.stats.txt" 