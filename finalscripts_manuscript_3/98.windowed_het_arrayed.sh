#!/usr/bin/env bash

#SBATCH -J BEAST # A single job name for the array
#SBATCH --ntasks-per-node=2 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 1-0:00:00 ### 15 seconds
#SBATCH --mem 90G
#SBATCH -o /scratch/rjp5nc/erroroutputs/beast.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/erroroutputs/beast.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab
#SBATCH --mail-type=END               # Send email at job completion
#SBATCH --mail-user=rjp5nc@virginia.edu    # Email address for notifications
#SBATCH --array=1-12

module load bcftools

VCF="/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/trimmed10bp_allsites_usobtusa.bgz.vcf.gz"
WINDOW=100000
RESULTDIR="/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/windowed_hets100000_all"
mkdir -p "$RESULTDIR"

# Get contig for this array index
contig=$(sed -n "${SLURM_ARRAY_TASK_ID}p" contigs.txt)
echo "Processing contig: $contig"

# Get sample names
bcftools query -l "$VCF" > "$RESULTDIR/samples.txt"
mapfile -t samples < "$RESULTDIR/samples.txt"
nsamples=${#samples[@]}

# Export sample names for awk
for i in "${!samples[@]}"; do
    export "samples$i=${samples[$i]}"
done


# Main loop
bcftools query -f '%CHROM\t%POS[\t%GT]\n' -r "$contig" "$VCF" | \
awk -v nsamples="$nsamples" -v contig="$contig" -v win_size="$WINDOW" '
{
    win = int($2 / win_size)
    for(i=3;i<=NF;i++){
        key = contig ":" win ":" i
        if($i ~ /^0[\/|]1$/ || $i ~ /^1[\/|]0$/){
            het[key]++
        }
        count[key]++
    }
}
END {
    for(k in count){
        split(k,a,":")
        win_start = a[2]*win_size
        win_end = win_start + win_size - 1
        sample_name = ENVIRON["samples" a[3]-3]
        het_prop = (count[k]>0) ? het[k]/count[k] : 0
        printf "%s\t%d\t%d\t%s\t%.5f\n", a[1], win_start, win_end, sample_name, het_prop
    }
}' | sort -k2,2n -k4,4 > "$RESULTDIR/${contig}_het_100kb.txt"

echo "Done: $contig"

#cat $RESULTDIR/*_het_10kb.txt > /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/only2_het_10kb_all.txt