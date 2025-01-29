
#!/usr/bin/env bash
#
#SBATCH -J toFastq # A single job name for the array
#SBATCH --ntasks-per-node=10 # ten core
#SBATCH -N 1 # on one node
#SBATCH -t 0-5:00  ### 48 hours
#SBATCH --mem 60G
#SBATCH -o /scratch/rjp5nc/basicstats.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/basicstats.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# sbatch ~/Genoseq_2023/basicStats.sh
### sacct -j 45333345

### Daphnia reference mitochondrial genome
### Scaffold_7757_HRSCAF_8726 6490423 to 6492626

module load samtools

samtools fastq /project/berglandlab/Robert/HMWDNAElvis3/m84128_250121_222443_s2.hifi_reads.bc2104.bam > /project/berglandlab/Robert/HMWDNAElvis3/fastq/Elvis3.fastq
