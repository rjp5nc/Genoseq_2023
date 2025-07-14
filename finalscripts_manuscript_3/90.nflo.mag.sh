#!/usr/bin/env bash
#SBATCH -J nfloMag    # Job name
#SBATCH --ntasks=1        # Single task per job
#SBATCH --cpus-per-task=10 # Number of CPU cores per task
#SBATCH -N 1              # Run on one node
#SBATCH -t 5-20:00        # 10 hours runtime
#SBATCH --mem=100G        # Memory per node
#SBATCH -o /scratch/rjp5nc/outputerrors/mag.%A_%a.out  # Standard output
#SBATCH -e /scratch/rjp5nc/outputerrors/mag.%A_%a.err  # Standard error
#SBATCH -p standard       # Partition
#SBATCH --account=berglandlab

#/home/rjp5nc/.nextflow/assets/nf-core/mag/

#curl -s https://get.nextflow.io | bash
#chmod +x nextflow
#export PATH=~/bin:$PATH
#nextflow -v


#https://www.youtube.com/watch?v=IiorfDHeoLo

module load java/21
module load bioconda/py3.10

/home/rjp5nc/bin/nextflow run /home/rjp5nc/.nextflow/assets/nf-core/mag/ \
-profile conda \
--skip_nanoplot \
--bin_min_size 5000 \
--bin_max_size 10000000 \
--input /scratch/rjp5nc/mag/samples.csv \
--outdir /scratch/rjp5nc/mag 

