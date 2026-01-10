#!/usr/bin/env bash
#
#SBATCH -J pulex2 # A single job name for the array
#SBATCH --cpus-per-task=10 # ten cores
#SBATCH -N 1 # on one node
#SBATCH -t 16:00:00 ### 8 hours
#SBATCH --mem 60G
#SBATCH -o /scratch/rjp5nc/outputerrors/download_%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/outputerrors/download_%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab
#SBATCH --array=1-54   # Adjust the range based on the number of folders


#  OUT="/scratch/rjp5nc/rawdata/mysamps_ids_dpulex_europe.txt"

# awk -F',' '
#   NR>1 {
#     srr=$3; sp=$5; cont=$9
#     gsub(/"/,"",srr); gsub(/"/,"",sp); gsub(/"/,"",cont)
#     if (sp=="Daphnia pulex" && cont=="Europe") print srr
#   }
# ' "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/Sample_ID_Species_merged_20251227.csv" | sort -u > "$OUT"

# echo "Wrote: $OUT"
# wc -l "$OUT"
# head "$OUT"

# SAMPLE_PARENT="/scratch/rjp5nc/rawdata/mysamples"
# IDFILE="/scratch/rjp5nc/rawdata/mysamps_ids_dpulex_europe.txt"

# tmp="${IDFILE}.tmp"
# : > "$tmp"

# while read -r id; do
#   id=$(echo "$id" | tr -d '\r' | xargs)
#   [[ -z "$id" ]] && continue

#   d=$(find "$SAMPLE_PARENT" -mindepth 1 -maxdepth 1 -type d -name "*${id}*" 2>/dev/null | head -n 1 || true)
#   if [[ -n "$d" ]]; then
#     echo "$id" >> "$tmp"
#   else
#     echo "DROP (no dir): $id" >&2
#   fi
# done < "$IDFILE"

# sort -u "$tmp" > "$IDFILE"
# rm -f "$tmp"

# echo "Cleaned IDFILE: $IDFILE"
# wc -l "$IDFILE"
# head "$IDFILE"


SAMPLE_PARENT="/scratch/rjp5nc/rawdata/mysamples"
IDFILE="/scratch/rjp5nc/rawdata/mysamps_ids_dpulex_europe.txt"

outroot="/scratch/rjp5nc/UK2022_2024/redone_mito/eupulex/mysamps"
DIRLIST="${outroot}/eu_dpulex_dirs.txt"
mkdir -p "$outroot"

# # build new dirlist safely
# tmp="${DIRLIST}.tmp"
# : > "$tmp"

# while read -r id; do
#   id=$(echo "$id" | tr -d '\r' | xargs)
#   [[ -z "$id" ]] && continue

#   d=$(find "$SAMPLE_PARENT" -mindepth 1 -maxdepth 1 -type d -name "*${id}*" 2>/dev/null | sort | head -n 1 || true)
#   if [[ -n "$d" ]]; then
#     printf "%s\n" "$d" >> "$tmp"
#   else
#     echo "WARN: no directory for ID=$id" >&2
#   fi
# done < "$IDFILE"

# sort -u "$tmp" > "$DIRLIST"
# rm -f "$tmp"

# echo "Wrote: $DIRLIST"
# wc -l "$DIRLIST"
# head "$DIRLIST"




module load cutadapt gcc/11.4.0 bwa/0.7.17 samtools/1.17 picard/2.27.5

SAMPLE_PARENT="/scratch/rjp5nc/rawdata/mysamples"
REF="/scratch/rjp5nc/Reference_genomes/mito_reference/"

outroot="/scratch/rjp5nc/UK2022_2024/redone_mito/eupulex/mysamps"
trimdir="${outroot}/trimmed_fastq"
bamdir="${outroot}/lane_bams"
mergedir="${outroot}/mergedbams"
sortdir="${outroot}/sortedbams"
repdir="${outroot}/sortedbamsreport"

mkdir -p "${trimdir}" "${bamdir}" "${mergedir}" "${sortdir}" "${repdir}"

# -----------------------
# 1) Build SRR list: European D. pulex (bash only)
# -----------------------
IDFILE="/scratch/rjp5nc/rawdata/mysamps_ids_dpulex_europe.txt"

# -----------------------
# 2) Map SRRs -> sample directories under SAMPLE_PARENT
#    Assumption: directory name contains the SRR somewhere.
#    (If not true on your system, tell me the naming pattern and I’ll adjust.)
# -----------------------
DIRLIST="/scratch/rjp5nc/UK2022_2024/redone_mito/eupulex/mysamps/eu_dpulex_dirs.txt"

# pick this dir
dir=$(awk -v n="${SLURM_ARRAY_TASK_ID}" 'NR==n{gsub(/\r/,""); print; exit}' "$DIRLIST")

if [[ -z "${dir}" ]]; then
  echo "Error: could not resolve sample dir for task ${SLURM_ARRAY_TASK_ID}"
  echo "DIRLIST: $DIRLIST"
  echo "DIRLIST lines: $(wc -l < "$DIRLIST")"
  exit 1
fi

if [[ ! -d "$dir" ]]; then
  echo "Error: dir not found on node: $dir"
  exit 1
fi

sample_name="$(basename "${dir}")"
echo "SAMPLE_DIR: ${dir}"
echo "SAMPLE:     ${sample_name}"
echo "REF:        ${REF}"

# -----------------------
# LOCATE LANE FASTQs INSIDE THIS SAMPLE DIR
# -----------------------
L6_1=$(ls "${dir}"/*L6_1.fq.gz 2>/dev/null | head -n 1 || true)
L6_2=$(ls "${dir}"/*L6_2.fq.gz 2>/dev/null | head -n 1 || true)
L3_1=$(ls "${dir}"/*L3_1.fq.gz 2>/dev/null | head -n 1 || true)
L3_2=$(ls "${dir}"/*L3_2.fq.gz 2>/dev/null | head -n 1 || true)
L4_1=$(ls "${dir}"/*L4_1.fq.gz 2>/dev/null | head -n 1 || true)
L4_2=$(ls "${dir}"/*L4_2.fq.gz 2>/dev/null | head -n 1 || true)

if [[ -z "${L6_1}" || -z "${L6_2}" || -z "${L3_1}" || -z "${L3_2}" || -z "${L4_1}" || -z "${L4_2}" ]]; then
  echo "Error: missing one or more lane FASTQs in ${dir}"
  echo "L6_1=${L6_1}"
  echo "L6_2=${L6_2}"
  echo "L3_1=${L3_1}"
  echo "L3_2=${L3_2}"
  echo "L4_1=${L4_1}"
  echo "L4_2=${L4_2}"
  ls -lh "${dir}" || true
  exit 1
fi

# -----------------------
# HELPER: trim + map one lane
# -----------------------
trim_and_map () {
  local lane="$1"
  local r1="$2"
  local r2="$3"

  local t1="${trimdir}/${sample_name}.${lane}.trimmed1.fq.gz"
  local t2="${trimdir}/${sample_name}.${lane}.trimmed2.fq.gz"
  local obam="${bamdir}/${sample_name}.${lane}.bam"

  echo "== Lane ${lane} =="
  echo "R1=${r1}"
  echo "R2=${r2}"

  cutadapt \
    -q 18 \
    --minimum-length 75 \
    -O 15 \
    -n 3 \
    --cores="${SLURM_CPUS_PER_TASK}" \
    -o "${t1}" \
    -p "${t2}" \
    "${r1}" "${r2}"

  bwa mem \
    -t "${SLURM_CPUS_PER_TASK}" \
    -R "@RG\tID:${sample_name}_${lane}\tSM:${sample_name}\tPL:illumina\tLB:lib1" \
    "${REF}" \
    "${t1}" "${t2}" | \
  samtools view -@ "${SLURM_CPUS_PER_TASK}" -Sbh -q 20 -F 0x100 - > "${obam}"
}

# -----------------------
# RUN LANES SEPARATELY
# -----------------------
trim_and_map "L6" "${L6_1}" "${L6_2}"
trim_and_map "L3" "${L3_1}" "${L3_2}"
trim_and_map "L4" "${L4_1}" "${L4_2}"

# -----------------------
# MERGE -> SORT -> DEDUP
# -----------------------
merged="${mergedir}/${sample_name}.merged.bam"
sorted="${sortdir}/${sample_name}.sorted.bam"

java -jar "$EBROOTPICARD/picard.jar" MergeSamFiles \
  -I "${bamdir}/${sample_name}.L3.bam" \
  -I "${bamdir}/${sample_name}.L4.bam" \
  -I "${bamdir}/${sample_name}.L6.bam" \
  -O "${merged}" \
  --SORT_ORDER unsorted

java -jar "$EBROOTPICARD/picard.jar" SortSam \
  -I "${merged}" \
  -O "${sorted}" \
  -SORT_ORDER coordinate

final_bam="/scratch/rjp5nc/UK2022_2024/redone_mito/eupulex/sortedbamsdedup_pulex/${sample_name}.sort.dedup.bam"
metrics="${repdir}/${sample_name}.mark_duplicates_report.txt"

java -jar "$EBROOTPICARD/picard.jar" MarkDuplicates \
  -REMOVE_DUPLICATES true \
  -I "${sorted}" \
  -O "${final_bam}" \
  -M "${metrics}" \
  -VALIDATION_STRINGENCY SILENT \
  -CREATE_INDEX true

echo "FINAL BAM: ${final_bam}"
echo "METRICS:   ${metrics}"
echo "Done: ${sample_name}"
