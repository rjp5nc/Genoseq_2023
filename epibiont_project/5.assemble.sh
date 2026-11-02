#!/usr/bin/env bash
#
#SBATCH -J megahit_v2
#SBATCH --ntasks-per-node=40
#SBATCH -N 1
#SBATCH -t 72:00:00
#SBATCH --mem 200G
#SBATCH -o /scratch/rjp5nc/err/megahit_v2.%A_%a.out
#SBATCH -e /scratch/rjp5nc/err/megahit_v2.%A_%a.err
#SBATCH -p standard
#SBATCH --account berglandlab_standard
#SBATCH --mail-type=END
#SBATCH --mail-user=rjp5nc@virginia.edu

set -euo pipefail
module load miniforge
export PATH=$HOME/MEGAHIT-1.2.9-Linux-x86_64-static/bin:$PATH

cd /scratch/rjp5nc/vorticella

megahit -1 nonDaphnia2_R1.fq.gz -2 nonDaphnia2_R2.fq.gz \
        -o megahit_out2 \
        -t 40 \
        -m 0.9 \
        --min-contig-len 500





#     #!/usr/bin/env bash
# set -euo pipefail
# module load seqkit
# module load blast

# cd /scratch/rjp5nc/vorticella
# asm=megahit_out2/final.contigs.fa
# out=analysis_v2
# mkdir -p "$out"

# count_telo () {
#   seqkit fx2tab "$1" | awk '{s=$NF; h=substr(s,1,30); t=substr(s,length(s)-29)}
#     h ~ /(CCCCAA){3}|(CCCCAAAA){2}/ || t ~ /(TTGGGG){3}|(TTTTGGGG){2}/ {c++}
#     END {print c+0, "of", NR, "contigs have telomere repeats at an end"}'
# }

# echo "===== Assembly stats ====="
# seqkit stats -a "$asm"

# echo; echo "===== Adapter check (should be near 0) ====="
# seqkit grep -s -i -p CTGTCTCTTATACACATCT -p AGATGTGTATAAGAGACAG "$asm" | grep -c ">" || true

# echo; echo "===== GC distribution ====="
# seqkit fx2tab -n -g -l "$asm" \
#   | awk '{gc=$NF; len=$(NF-1); b=int(gc/5)*5; n[b]++; s[b]+=len}
#          END {for (k in n) printf "%d-%d%%\t%d contigs\t%d bp\n", k, k+5, n[k], s[k]}' \
#   | sort -n

# echo; echo "===== Low-GC (<35%) and mid-GC (40-50%) contigs >= 1 kb ====="
# seqkit fx2tab -n -g -l "$asm" | awk '$NF < 35 && $(NF-1) >= 1000 {print $1}' > "$out"/lowGC_ids.txt
# seqkit fx2tab -n -g -l "$asm" | awk '$NF >= 40 && $NF < 50 && $(NF-1) >= 1000 {print $1}' > "$out"/midGC_ids.txt
# seqkit grep -f "$out"/lowGC_ids.txt "$asm" > "$out"/lowGC_1kb.fa
# seqkit grep -f "$out"/midGC_ids.txt "$asm" > "$out"/midGC_1kb.fa
# seqkit stats "$out"/lowGC_1kb.fa "$out"/midGC_1kb.fa

# echo; echo "===== Telomere counts ====="
# echo -n "low GC: "; count_telo "$out"/lowGC_1kb.fa
# echo -n "mid GC: "; count_telo "$out"/midGC_1kb.fa

# echo; echo "===== Telomere type in low-GC set ====="
# seqkit fx2tab "$out"/lowGC_1kb.fa | awk '{s=$NF; h=substr(s,1,30); t=substr(s,length(s)-29)}
#   {if (h ~ /(CCCCAAAA){2}/ || t ~ /(TTTTGGGG){2}/) a++;
#    else if (h ~ /(CCCCAA){3}/ || t ~ /(TTGGGG){3}/) b++}
#   END {print "C4A4:", a+0, "  C4A2:", b+0}'

# echo; echo "===== Complete chromosomes (telomeres at both ends) ====="
# seqkit fx2tab "$out"/lowGC_1kb.fa | awk '{s=$NF; h=substr(s,1,30); t=substr(s,length(s)-29)}
#   (h ~ /(CCCCAA){3}|(TTGGGG){3}/) && (t ~ /(TTGGGG){3}|(CCCCAA){3}/) {print $1, length(s)}' \
#   > "$out"/complete_chromosomes.txt
# echo "count: $(wc -l < "$out"/complete_chromosomes.txt)"
# if [ -s "$out"/complete_chromosomes.txt ]; then
#   awk '{print $2}' "$out"/complete_chromosomes.txt | sort -n \
#     | awk '{a[NR]=$1} END {print "min:", a[1], " median:", a[int(NR/2)+1], " max:", a[NR]}'
#   seqkit grep -f <(cut -d' ' -f1 "$out"/complete_chromosomes.txt) "$asm" > "$out"/complete_chr.fa
# fi

# echo; echo "===== rRNA BLAST ====="
# makeblastdb -in "$asm" -dbtype nucl -out "$out"/asm_db > /dev/null
# blastn -query RP03-Premixed.fa -db "$out"/asm_db \
#   -outfmt "6 qseqid sseqid pident length qlen slen evalue bitscore" \
#   -evalue 1e-20 -max_target_seqs 10 > "$out"/RP03_hits.tsv
# cat "$out"/RP03_hits.tsv

# top=$(sort -k8,8nr "$out"/RP03_hits.tsv | head -1 | cut -f2)
# echo; echo "Top rDNA contig: $top"
# seqkit grep -p "$top" "$asm" > "$out"/vort_rDNA.fa
# grep ">" "$out"/vort_rDNA.fa
# seqkit fx2tab -n -g -l "$out"/vort_rDNA.fa
# echo "Contig ends:"
# seqkit fx2tab "$out"/vort_rDNA.fa | awk '{s=$NF; print substr(s,1,40); print substr(s,length(s)-39)}'