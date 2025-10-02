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
#SBATCH --array=1-12

#cat /scratch/rjp5nc/err/pixy.4130404_8

cd /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/

module load bcftools
#bcftools index /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/trimmed10bp_allsites_usobtusa.bgz.vcf.gz

VCF=trimmed10bp_allsites_usobtusa.bgz.vcf.gz

 bcftools query -l /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/trimmed10bp_allsites_usobtusa.bgz.vcf.gz > samples.txt
nsamples=$(bcftools query -l "$VCF" | wc -l)
echo "Samples: $nsamples"

samples_str=$(bcftools query -l "$VCF" | paste -sd "," -)
echo "Samples: $samples_str"

RESULTDIR=dp_results_100k_unfiltered
mkdir -p "$RESULTDIR"

# Get contig for this array task
contig=$(sed -n "${SLURM_ARRAY_TASK_ID}p" contigs_clean.txt)
if [[ -z "$contig" ]]; then
    echo "ERROR: contig is empty for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"
    exit 1
fi
echo "Processing contig: $contig"

bcftools query -f '%CHROM\t%POS[\t%DP]\n' -r "$contig" "$VCF" \
| awk -v nsamples="$nsamples" -v contig="$contig" -v SAMPLES="$samples_str" '
BEGIN { split(SAMPLES, sample_arr, ",") }
{
    win = int($2/100000)
    for(i=1;i<=nsamples;i++){
        dp = $(i+2)
        if(dp ~ /^[0-9]+$/){ 
            sum[i,win] += dp
            count[i,win]++
        }
        # mark the window even if DP missing
        seen[win]=1
    }
}
END {
    for(win in seen){
        win_start = win*100000
        win_end = win_start + 99999
        for(i=1;i<=nsamples;i++){
            sample_name = sample_arr[i]
            avg = (count[i,win]>0) ? sum[i,win]/count[i,win] : 0
            print contig, win_start, win_end, sample_name, avg
        }
    }
}' \
| sort -k2,2n -k4,4 \
> "$RESULTDIR/${contig}.avgdepth.long.sorted_100000.txt"

# cat *avgdepth.long.sorted_10000.txt > ../avgdepth_all_unfilt_10000.txt

# awk '
# {
#     key = $1"\t"$2"\t"$3
#     sum[key] += $5
#     count[key] += 1
# }
# END {
#     OFS="\t"
#     print "contig","window_start","window_end","avg_depth"
#     for(k in sum){
#         avg = sum[k]/count[k]
#         print k, avg
#     }
# }
# ' avgdepth_all_unfilt_10000.txt > window_avg_avgdepth_all_unfilt_10000.txt