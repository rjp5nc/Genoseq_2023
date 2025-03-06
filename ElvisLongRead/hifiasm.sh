#!/usr/bin/env bash
#
#SBATCH -J hifiasm # A single job name for the array
#SBATCH --ntasks-per-node=40 # multi core
#SBATCH -N 1 # on one node
#SBATCH -t 6-00:00 # 6 days
#SBATCH --mem 100G
#SBATCH -o /scratch/rjp5nc/Canu_error/hifiasm%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/Canu_error/hifiasm%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab


#conda create -n hifiasm -y
#conda activate hifiasm
#conda install -c bioconda hifiasm

conda activate hifiasm

wd="/scratch/rjp5nc/HMW/HMWDNAElvis3"
cd ${wd}

# Set variables
HIFI_READS="/scratch/rjp5nc/HMW/HMWDNAElvis3/m84128_250121_222443_s2.hifi_reads.bc2104.fq.gz"  # Change this to your HiFi reads file
THREADS=40  # Adjust based on your system
OUTPUT_PREFIX="assembly"

# Run HiFiasm
hifiasm -o ${OUTPUT_PREFIX} -t ${THREADS} ${HIFI_READS}

# Extract the primary haploid assembly
awk '/^>/ {if (seen[$1]++) next} {print}' ${OUTPUT_PREFIX}.p_ctg.fa > haploid_assembly.fa

# Done
echo "Hifiasm assembly complete! Primary contigs saved as haploid_assembly.fa"
