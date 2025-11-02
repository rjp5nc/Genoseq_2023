#!/usr/bin/env bash

#SBATCH -J BEAST # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 6-0:00:00 ### 15 seconds
#SBATCH --mem 60G
#SBATCH -o /scratch/rjp5nc/erroroutputs/beast.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/erroroutputs/beast.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab
#SBATCH --mail-type=END               # Send email at job completion
#SBATCH --mail-user=rjp5nc@virginia.edu    # Email address for notifications
#SBATCH --array=1-79


# Chromosome
uniqueid=$(sed -n "${SLURM_ARRAY_TASK_ID}p" /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/heterozygotegilmer/onlygilmer.txt)
BAM="/scratch/rjp5nc/UK2022_2024/final_bam_rg2/${uniqueid}finalmap_RG.bam"
OUTDIR="/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/depths"
mkdir -p $OUTDIR

/scratch/rjp5nc/mosdepth  \
  --by 100000 \
  --threads 4 \
  ${OUTDIR}/${uniqueid} \
  ${BAM}


  MERGED="all_samples_100kb_depths.tsv"

# header
echo -e "sample\tchrom\tstart\tend\tdepth" > /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/$MERGED

for f in $OUTDIR/*.regions.bed.gz
do
    sample=$(basename $f .regions.bed.gz)
    zcat $f | awk -v s=$sample '{print s"\t"$1"\t"$2"\t"$3"\t"$4}' >> /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/$MERGED
done
