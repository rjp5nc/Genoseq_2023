#!/usr/bin/env bash
#
#SBATCH -J busco_vort
#SBATCH --ntasks-per-node=16
#SBATCH -N 1
#SBATCH -t 10:00:00
#SBATCH --mem 100G
#SBATCH -o /scratch/rjp5nc/err/busco_vort.%A_%a.out
#SBATCH -e /scratch/rjp5nc/err/busco_vort.%A_%a.err
#SBATCH -p standard
#SBATCH --account berglandlab_standard
#SBATCH --mail-type=END
#SBATCH --mail-user=rjp5nc@virginia.edu

cd /scratch/rjp5nc/vorticella
module load bcftools

ls persample/Rockpool4_*.bam > persample/bams.txt
awk '{print $1"\t1\t10000000"}' vort_bin/vort_bin2_ids.txt > vort_bin/vort_bin2_regions.tsv

bcftools mpileup -f megahit_out2/final.contigs.fa -R vort_bin/vort_bin2_regions.tsv \
  -b persample/bams.txt -a AD,DP -q 20 -Q 20 --threads 16 \
  | bcftools call -mv --ploidy 1 --threads 16 -Oz -o vort_bin/vort_persample_snps.vcf.gz
bcftools index vort_bin/vort_persample_snps.vcf.gz

# Number of SNPs and SNP density
n=$(bcftools view -H -i 'QUAL>=30' vort_bin/vort_persample_snps.vcf.gz | wc -l)
echo "SNPs: $n  (per 10 kb: $(awk -v n=$n 'BEGIN{printf "%.2f", n/1000}'))"