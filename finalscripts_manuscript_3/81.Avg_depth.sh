#!/usr/bin/env bash
#
#SBATCH -J pixy # A single job name for the array
#SBATCH --ntasks-per-node=10 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-92:00  ### 48 hours
#SBATCH --mem 30G
#SBATCH -o /scratch/rjp5nc/err/pixy.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/err/pixy.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab


cd /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/

module load bcftools

VCF=/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/trimmed10bp_allsites_Repeatmasked_usobtusa.filtered_bgz2.vcf.gz

# pick the contig for this task
contig=$(sed -n "${SLURM_ARRAY_TASK_ID}p" contigs.txt)

# run bcftools query
bcftools query -r $contig -f '%CHROM\t%POS[\t%DP]\n' $VCF \
| awk '{
    win = int($2/10000);
    key = $1":"win;
    for (i=3; i<=NF; i++) {
      sum[key,i] += $i;
      count[key,i]++;
    }
  }
  END {
    for (k in sum) {
      split(k,a,SUBSEP);
      chromwin=a[1];
      sample=a[2]-2;
      avg = sum[k]/count[k];
      print chromwin, sample, avg;
    }
  }' > results/${contig}.avgdepth.txt


  #cat results/*.avgdepth.txt > all_samples.avgdepth.txt
