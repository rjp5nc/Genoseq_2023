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

#cat /scratch/rjp5nc/err/pixy.4130252_4

cd /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/

module load bcftools

VCF=/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/trimmed10bp_allsites_Repeatmasked_usobtusa.filtered_bgz.vcf.gz

bcftools query -l /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/trimmed10bp_allsites_Repeatmasked_usobtusa.filtered_bgz.vcf.gz \
> samples.txt

bcftools index -s /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/trimmed10bp_allsites_Repeatmasked_usobtusa.filtered_bgz.vcf.gz \
| cut -f1 > contigs.txt

sed -i 's/\r//g; s/ //g' contigs.txt

cut -f1 contigs.txt > contigs_only.txt
grep -v "^$" contigs_only.txt > contigs_clean.txt

VCF=/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/trimmed10bp_allsites_Repeatmasked_usobtusa.filtered_bgz.vcf.gz
RESULTDIR=results
mkdir -p "$RESULTDIR"

# get contig for this array task
contig=$(sed -n "${SLURM_ARRAY_TASK_ID}p" contigs.txt)

# get sample names into an array
mapfile -t samples < samples.txt
nsamples=${#samples[@]}
echo $nsamples  # optional check

# export sample names to ENV for AWK
for i in "${!samples[@]}"; do
    export "samples$i=${samples[$i]}"
done

# bcftools query -f '%CHROM\t%POS[\t%DP]\n' -r "$contig" "$VCF" \
# | awk -v nsamples="$nsamples" -v contig="$contig" '
# {
#     win = int($2/10000)
#     for(i=3;i<=NF;i++){
#         key=contig":"win":"i
#         sum[key] += ($i!=".") ? $i : 0
#         count[key] += ($i!=".")
#     }
# }
# END {
#     for(k in sum){
#         split(k,a,":")
#         win_start = a[2]*10000
#         win_end = win_start + 9999
#         sample_name = ENVIRON["samples" a[3]-3]
#         avg = (count[k]>0) ? sum[k]/count[k] : 0
#         printf "%s\t%d\t%d\t%s\t%.2f\n", a[1], win_start, win_end, sample_name, avg
#     }
# }' \
# | sort -k2,2n -k4,4 \
# > "$RESULTDIR/${contig}.avgdepth.long.sorted.txt"

#cat *.avgdepth.long.sorted.txt > ../avg.depth.all12.txt


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
# ' avg.depth.all12.txt > window_avg_depth.txt





# get contig for this array task
contig=$(sed -n "${SLURM_ARRAY_TASK_ID}p" contigs_clean.txt)

# get sample names into an array
mapfile -t samples < samples.txt
nsamples=${#samples[@]}

# export sample names to ENV for AWK
for i in "${!samples[@]}"; do
    export "samples$i=${samples[$i]}"
done


# main depth extraction + window averaging
bcftools query -f '%CHROM\t%POS[\t%DP]\n' -r "$contig" "$VCF" \
| awk -v nsamples="$nsamples" -v contig="$contig" '
{
    win = int($2/100000)   # 100kb windows
    for(i=3;i<=NF;i++){
        key=contig":"win":"i
        sum[key]   += ($i!=".") ? $i : 0
        count[key] += ($i!=".")
    }
}
END {
    for(k in sum){
        split(k,a,":")
        win_start = a[2]*100000
        win_end   = win_start + 99999
        sample_name = ENVIRON["samples" a[3]-3]
        avg = (count[k] > 0) ? sum[k] / count[k] : 0
        # always print
        printf "%s\t%d\t%d\t%s\t%.2f\n", a[1], win_start, win_end, sample_name, avg
    }
}' \
| sort -k2,2n -k4,4 \
> "$RESULTDIR/${contig}.avgdepth.long.sorted_100000.txt"




#cat *.avgdepth.long.sorted_100000.txt > ../avg.depth.all12_100000.txt

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
# ' avg.depth.all12_100000.txt > window_avg_depth_100000.txt