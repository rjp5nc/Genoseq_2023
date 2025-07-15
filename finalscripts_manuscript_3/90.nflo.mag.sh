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

#module load java/21
#module load bioconda/py3.10


#conda create -n nfcore-mag-env python=3.10 pandas=1.4.3 plotly=5.13 nanoplot=1.41.6 -c conda-forge -c bioconda -y
#conda install -n nfcore-mag-env multiqc samtools nextflow -c bioconda -c conda-forge
#conda install -n nfcore-mag-env nanoplot=1.40.0 -c bioconda -c conda-forge

source ~/miniconda3/etc/profile.d/conda.sh
conda activate nfcore-mag-env


cd /scratch/rjp5nc/mag


#NanoPlot \
#  -p raw_ --title Elvis03_raw -c darkblue \
#  -t 2 \
#  --fastq /project/berglandlab/Robert/HMWDNAElvis3/RawDatafromCDGenomics/fastq/m84128_250121_222443_s2.hifi_reads.bc2104.fq.gz

/home/rjp5nc/bin/nextflow run /home/rjp5nc/.nextflow/assets/nf-core/mag/ \
  -profile conda \
  --input /scratch/rjp5nc/mag/samples.csv \
  --outdir /scratch/rjp5nc/mag \
  --validate_params false


#--bin_min_size 5 \
#--bin_max_size 100000000 \


