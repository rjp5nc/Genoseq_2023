#!/usr/bin/env bash
#
#SBATCH -J vort_mito_check
#SBATCH --ntasks-per-node=16
#SBATCH -N 1
#SBATCH -t 12:00:00
#SBATCH --mem 32G
#SBATCH -o /scratch/rjp5nc/err/vort_mito_check.%A_%a.out
#SBATCH -e /scratch/rjp5nc/err/vort_mito_check.%A_%a.err
#SBATCH -p standard
#SBATCH --account berglandlab_standard
#SBATCH --mail-type=END
#SBATCH --mail-user=rjp5nc@virginia.edu

set -euo pipefail
module purge
module load samtools bcftools

cd /scratch/rjp5nc/vorticella

# 1. Read-pair links between contigs (pairs with mates on different contigs, MAPQ>=20, >=3 pairs)
if [ ! -s vort_bin/contig_links.tsv ]; then
  samtools view -@ 8 -q 20 -F 0x904 -f 1 coverage_v2/contigs_v2.sorted.bam \
    | awk '$7!="=" && $7!="*" {print $3"\t"$7}' \
    | awk '{k=($1<$2)?$1"\t"$2:$2"\t"$1; n[k]++} END {for (x in n) if (n[x]>=3) print x"\t"n[x]}' \
    > vort_bin/contig_links.tsv
fi
echo "Contig links: $(wc -l < vort_bin/contig_links.tsv)"

# 2. Contigs in the mito depth/GC window linked to the known mito contigs
awk -F'\t' 'FILENAME==ARGV[1] {seed[$1]=1; next}
            FILENAME==ARGV[2] {win[$1]=1; next}
  ($1 in seed) && !($2 in seed) && ($2 in win) {print $2"\t"$1"\t"$3}
  ($2 in seed) && !($1 in seed) && ($1 in win) {print $1"\t"$2"\t"$3}' \
  vort_bin/vort_mito_ids.txt vort_bin/mito_window_ids.txt vort_bin/contig_links.tsv \
  | sort -k3,3nr > vort_bin/mito_linked.tsv

echo "== Linked mito candidates (contig, linked_to, n_pairs, length, depth, gc) =="
awk -F'\t' 'NR==FNR {l[$1]=$2"\t"$3"\t"$4; next} {print $0"\t"l[$1]}' \
  coverage_v2/contig_length_depth_gc.tsv vort_bin/mito_linked.tsv

# 3. Mito variants across the 36 samples (known mito contigs)
awk '{print $1"\t1\t10000000"}' vort_bin/vort_mito_ids.txt > vort_bin/vort_mito_regions.tsv
bcftools mpileup -f megahit_out2/final.contigs.fa -R vort_bin/vort_mito_regions.tsv \
  -b persample/bams.txt -a AD,DP -d 10000 -q 20 -Q 20 --threads 16 \
  | bcftools call -mv --ploidy 1 --threads 16 -Oz -o vort_bin/vort_mito_snps.vcf.gz
bcftools index -f vort_bin/vort_mito_snps.vcf.gz

echo "Mito SNPs (QUAL>=30): $(bcftools view -H -i 'QUAL>=30' -v snps vort_bin/vort_mito_snps.vcf.gz | wc -l)"

bcftools view -i 'QUAL>=30' -v snps -m2 -M2 vort_bin/vort_mito_snps.vcf.gz \
  | bcftools query -H -f '%CHROM\t%POS[\t%AD]\n' \
  | awk -F'\t' 'NR==1 {printf "contig\tpos"; for (i=3;i<=NF;i++) {n=$i; sub(/.*Rockpool4_/,"",n); sub(/\.bam.*/,"",n); printf "\t%s", n}; print ""; next}
    { printf "%s\t%s", $1, $2
      for (i=3;i<=NF;i++) { split($i,a,","); d=a[1]+a[2]; printf "\t%s", (d>=20) ? sprintf("%.2f", a[2]/d) : "NA" }
      print "" }' > vort_bin/vort_mito_af.tsv

echo "== Pure vs mixed per sample-site (pure: AF<0.1 or >0.9) =="
awk -F'\t' 'NR>1 {for (i=3;i<=NF;i++) if ($i!="NA") {
    c = ($i<0.1 || $i>0.9) ? "pure" : "mixed"; n[c]++ }}
  END {for (k in n) print k"\t"n[k]}' vort_bin/vort_mito_af.tsv