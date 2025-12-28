#!/usr/bin/env bash
#SBATCH -J DownloadMap    # Job name
#SBATCH --ntasks=1        # Single task per job
#SBATCH --cpus-per-task=10 # Number of CPU cores per task
#SBATCH -N 1              # Run on one node
#SBATCH -t 1-20:00        # 10 hours runtime
#SBATCH --mem=50G        # Memory per node
#SBATCH -o /scratch/rjp5nc/outputerrors/down.%A_%a.out  # Standard output
#SBATCH -e /scratch/rjp5nc/outputerrors/down.%A_%a.err  # Standard error
#SBATCH -p standard       # Partition
#SBATCH --account=berglandlab


cd /scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/
module load bcftools samtools

reference="/scratch/rjp5nc/Reference_genomes/mito_reference/usdobtusa_mito.fasta"
outdir="/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito"
vcf="${outdir}/usdobtusa_mito_allsites_all.diploid.vcf.gz"

mkdir -p "$outdir"
[ -f "${reference}.fai" ] || samtools faidx "$reference"

# Build BAM list (skip CSV header)
# awk -F, 'NR>1 {gsub(/[^a-zA-Z0-9_]/,"",$6); if ($6=="usdobtusa_mito") print $1 "finalmap_RG.bam"}' \
#   "$meta" > usdobtusa_bams.txt


find /scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/final_mitobam_rg3/ \
  -name "*finalmap_RG.bam" \
| xargs -n 1 basename \
| sort \
> /scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/usdobtusa_bams.txt



# Multi-sample call in one go, haploid, keep all covered sites
cd /scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/final_mitobam_rg3/
bcftools mpileup \
  -f "$reference" \
  -q 20 -Q 20 \
  -a FORMAT/DP,FORMAT/AD \
  -Ou -b /scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/usdobtusa_bams.txt \
| bcftools call -m -A --ploidy 2 -Oz -o "$vcf"








bcftools index -t "$vcf"
echo "Wrote haploid all-sites multi-sample VCF -> $vcf"
