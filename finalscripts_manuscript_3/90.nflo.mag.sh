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

#rm /scratch/rjp5nc/mag/samplesheet.csv
#cat <<EOF > /scratch/rjp5nc/mag/samplesheet.csv
#sample,group,long_reads
#Elvis3_HiFi,test_group,/project/berglandlab/Robert/HMWDNAElvis3/RawDatafromCDGenomics/fastq/m84128_250121_222443_s2.hifi_reads.bc2104.fq.gz
#EOF



/home/rjp5nc/bin/nextflow run nf-core/mag \
  -profile singularity \
  --input '/project/berglandlab/Robert/HMWDNAElvis3/RawDatafromCDGenomics/fastq/m84128_250121_222443_s2.hifi_reads.bc2104.fq.gz' \
  --outdir /scratch/rjp5nc/mag 