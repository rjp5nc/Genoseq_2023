#!/usr/bin/env bash
#
#SBATCH -J move # A single job name for the array
#SBATCH --ntasks-per-node=60 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0:24:00  ### 48 hours
#SBATCH --mem 60G
#SBATCH -o /scratch/rjp5nc/err/move.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/err/move.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab
#SBATCH --mail-type=END               # Send email at job completion
#SBATCH --mail-user=rjp5nc@virginia.edu    # Email address for notifications

set -euo pipefail

module load bowtie2
module load samtools

threads=60
wd=/scratch/rjp5nc/vorticella
ref=/scratch/rjp5nc/Reference_genomes/post_kraken/US_obtusa_onlydaps.fa
idxdir=$wd/daphnia_idx
idx=$idxdir/US_obtusa

mkdir -p "$idxdir"
cd "$wd"

# Build Daphnia index (skipped if it already exists)
if [ ! -e "${idx}.1.bt2" ] && [ ! -e "${idx}.1.bt2l" ]; then
  bowtie2-build --threads "$threads" "$ref" "$idx"
fi

# Map, keep pairs where both reads are unmapped, sort by name
bowtie2 -x "$idx" -1 trim_R1.fq.gz -2 trim_R2.fq.gz -p "$threads" \
        --very-sensitive 2> host_mapping.log \
  | samtools view -@ 4 -b -f 12 -F 256 - \
  | samtools sort -n -@ 4 -m 2G -T "$wd/sorttmp" -o nonDaphnia.namesorted.bam -

# Convert to paired FASTQ
samtools fastq -@ 4 -c 6 \
  -1 nonDaphnia_R1.fq.gz -2 nonDaphnia_R2.fq.gz \
  -0 /dev/null -s /dev/null -n \
  nonDaphnia.namesorted.bam

samtools flagstat nonDaphnia.namesorted.bam > nonDaphnia.flagstat.txt
cat host_mapping.log