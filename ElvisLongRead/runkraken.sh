#!/usr/bin/env bash
#
#SBATCH -J run_kraken # A single job name for the array
#SBATCH --ntasks-per-node=30 # multi core
#SBATCH -N 1 # on one node
#SBATCH -t 1-03:00 # 3 hours
#SBATCH --mem 150G
#SBATCH -o /scratch/rjp5nc/Canu_error/kraken%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/Canu_error/kraken%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab


module load kraken2

kraken2 --db /scratch/rjp5nc/krakenDB/nt --threads 16 --report /scratch/rjp5nc/krakenDB/report.txt --output /scratch/rjp5nc/krakenDB/kraken_output.txt /scratch/rjp5nc/HMW/HMWDNAElvis3/hifiasm_out/assembly.hap2.p_ctg.fa
