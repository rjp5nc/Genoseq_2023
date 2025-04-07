#!/usr/bin/env bash
#
#SBATCH -J liftover # A single job name for the array
#SBATCH --ntasks-per-node=10 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 2-10:00 # 10 hours
#SBATCH --mem 100G
#SBATCH -o /scratch/rjp5nc/erroroutputs/lift.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/erroroutputs/lift.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab


#conda create -n plotsr
conda deactivate
conda activate plotsr
#conda install numpy=1.21.2 pandas=1.2.4 matplotlib=3.3.4 setuptools
#conda install syri

#syri -c A_B.bam -r A.fa -q B.fa -F B --prefix A_B &

REF="eu_pulex_totalHiCwithallbestgapclosed.clean.fa"
QUERY="assembly.hap2_onlydaps.fasta"

cd /scratch/rjp5nc/Reference_genomes/liftover/

syri -c /scratch/rjp5nc/Reference_genomes/liftover/ribbon/euobtusa_v_eupulex.bam \
 -r /scratch/rjp5nc/Reference_genomes/post_kraken/${REF} \
 -q /scratch/rjp5nc/Reference_genomes/post_kraken/${QUERY} \
 -F B --prefix \
 euonbtusa_v_eupulex &

