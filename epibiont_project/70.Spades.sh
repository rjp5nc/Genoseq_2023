#!/usr/bin/env bash
#
#SBATCH -J vort_reassemble
#SBATCH --ntasks-per-node=40
#SBATCH -N 1
#SBATCH -t 48:00:00
#SBATCH --mem 120G
#SBATCH -o /scratch/rjp5nc/err/vort_reassemble.%A_%a.out
#SBATCH -e /scratch/rjp5nc/err/vort_reassemble.%A_%a.err
#SBATCH -p standard
#SBATCH --account berglandlab_standard
#SBATCH --mail-type=END
#SBATCH --mail-user=rjp5nc@virginia.edu

set -euo pipefail
module purge
module load gcc/14.2.0          # <- prerequisite from 'module spider spades' (adjust)
module load spades              # <- exact name/version from 'module spider spades' (adjust)
module load samtools seqkit

wd=/scratch/rjp5nc/vorticella
out=$wd/vort_reasm
mkdir -p "$out"
cd "$wd"

# 1. Contigs to recruit reads from: Vorticella window, >=500 bp
#    GC < 38%, depth 14-250x, not Kraken-bacterial/archaeal/viral/human
awk -F'\t' '
  FILENAME==ARGV[1] { kd[$1]=$2; kp[$1]=$3; next }
  FNR==1 { next }
  { id=$1; len=$2; dep=$3; gc=$4
    if (len < 500 || gc >= 38 || dep < 14 || dep > 250) next
    if (kd[id]=="Bacteria" || kd[id]=="Archaea" || kd[id]=="Viruses") next
    if (kd[id]=="Eukaryota" && kp[id]=="Chordata") next
    print id"\t0\t"len }' vort_bin/kraken_class.tsv coverage_v2/contig_length_depth_gc.tsv \
  > "$out"/recruit_contigs.bed
echo "Recruit contigs: $(wc -l < "$out"/recruit_contigs.bed)"

# 2. Names of reads with at least one mate on those contigs (MAPQ >= 10)
samtools view -@ 16 -q 10 -F 0x904 -M -L "$out"/recruit_contigs.bed \
  coverage_v2/contigs_v2.sorted.bam \
  | cut -f1 | sort -u -S 8G --parallel=16 > "$out"/read_names.txt
echo "Read names: $(wc -l < "$out"/read_names.txt)"

# 3. Pull both mates from the clean reads
seqkit grep -j 16 -f "$out"/read_names.txt nonDaphnia2_R1.fq.gz -o "$out"/vort_R1.fq.gz
seqkit grep -j 16 -f "$out"/read_names.txt nonDaphnia2_R2.fq.gz -o "$out"/vort_R2.fq.gz
seqkit stats "$out"/vort_R1.fq.gz "$out"/vort_R2.fq.gz

# 4. Assemble (single-genome mode)
spades.py -1 "$out"/vort_R1.fq.gz -2 "$out"/vort_R2.fq.gz \
  -o "$out"/spades -t 40 -m 120 -k 21,33,55,77,99,127
echo "SPAdes done: $(date)"

# 5. QC on scaffolds >= 1 kb
seqkit seq -m 1000 "$out"/spades/scaffolds.fasta > "$out"/vort_reasm_1kb.fa
echo "== Old bin vs reassembly =="
seqkit stats -a vort_bin/vort_bin2.fa "$out"/vort_reasm_1kb.fa

echo "== Telomeres in reassembly =="
seqkit fx2tab "$out"/vort_reasm_1kb.fa | awk '{s=$NF; h=substr(s,1,30); t=substr(s,length(s)-29)}
  {e1=(h ~ /(CCCCAA){3}|(TTGGGG){3}/); e2=(t ~ /(TTGGGG){3}|(CCCCAA){3}/)}
  e1||e2 {one++} e1&&e2 {both++}
  END {print "Telomere at >=1 end:", one+0, "  Both ends (complete chromosomes):", both+0}'

# 6. BUSCO (alveolata, ciliate code table 6)
module purge
module load apptainer/1.4.5
module load busco/6.1.0
apptainer exec -B /scratch $CONTAINERDIR/busco-6.1.0.sif busco \
  -i "$out"/vort_reasm_1kb.fa -m genome -l alveolata_odb12 -c 40 \
  --offline --download_path /scratch/rjp5nc/buscoDB \
  --metaeuk \
  --metaeuk_parameters="--translation-table=6" \
  --metaeuk_rerun_parameters="--translation-table=6" \
  --out_path "$out" -o busco_alveolata_tt6 -f
echo "== BUSCO =="
grep -h "C:" "$out"/busco_alveolata_tt6/short_summary*.txt