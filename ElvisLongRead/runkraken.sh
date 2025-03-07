#!/usr/bin/env bash
#
#SBATCH -J run_kraken # A single job name for the array
#SBATCH --ntasks-per-node=10 # multi core
#SBATCH -N 1 # on one node
#SBATCH -t 1-03:00 # 3 hours
#SBATCH --mem 40G
#SBATCH -o /scratch/rjp5nc/Canu_error/kraken%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/Canu_error/kraken%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab


module load kraken2

kraken2 --db /scratch/rjp5nc/krakenDB/nt --threads 16 --report report.txt --output kraken_output.txt /scratch/rjp5nc/HMW/HMWDNAElvis3/hifiasm_out/assembly.hap2.p_ctg.fa
