#!/usr/bin/env bash
#
#SBATCH -J trim_map_v2
#SBATCH --ntasks-per-node=40
#SBATCH -N 1
#SBATCH -t 48:00:00
#SBATCH --mem 60G
#SBATCH -o /scratch/rjp5nc/err/trim_map_v2.%A_%a.out
#SBATCH -e /scratch/rjp5nc/err/trim_map_v2.%A_%a.err
#SBATCH -p standard
#SBATCH --account berglandlab_standard
#SBATCH --mail-type=END
#SBATCH --mail-user=rjp5nc@virginia.edu

set -euo pipefail

module load trimmomatic/0.40
module load bowtie2
module load samtools

threads=40
wd=/scratch/rjp5nc/vorticella
idx=$wd/daphnia_idx/US_obtusa
jar=$(ls "$EBROOTTRIMMOMATIC"/*.jar | head -1)

cd "$wd"

# Combine TruSeq + Nextera adapters
cat "$EBROOTTRIMMOMATIC"/adapters/TruSeq3-PE-2.fa \
    "$EBROOTTRIMMOMATIC"/adapters/NexteraPE-PE.fa > adapters_all.fa

# 1. Trim
java -Xmx24G -jar "$jar" PE -threads 16 -phred33 \
  Rockpool4_D1-F12_R1.fq.gz Rockpool4_D1-F12_R2.fq.gz \
  trim2_R1.fq.gz unpaired2_R1.fq.gz \
  trim2_R2.fq.gz unpaired2_R2.fq.gz \
  ILLUMINACLIP:adapters_all.fa:2:30:10:2:True \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50 \
  2> trimmomatic2.log

# 2. Map to Daphnia, keep pairs where both reads are unmapped
bowtie2 -x "$idx" -1 trim2_R1.fq.gz -2 trim2_R2.fq.gz -p "$threads" \
        --very-sensitive 2> host_mapping2.log \
  | samtools view -@ 4 -b -f 12 -F 256 - \
  | samtools sort -n -@ 4 -m 2G -T "$wd/sorttmp2" -o nonDaphnia2.namesorted.bam -

# 3. Convert to paired FASTQ
samtools fastq -@ 4 -c 6 \
  -1 nonDaphnia2_R1.fq.gz -2 nonDaphnia2_R2.fq.gz \
  -0 /dev/null -s /dev/null -n \
  nonDaphnia2.namesorted.bam

samtools flagstat nonDaphnia2.namesorted.bam > nonDaphnia2.flagstat.txt
cat trimmomatic2.log host_mapping2.log