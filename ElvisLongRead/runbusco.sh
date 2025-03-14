#!/usr/bin/env bash
#
#SBATCH -J Busco # A single job name for the array
#SBATCH --ntasks-per-node=10 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-05:00 # 5 hours
#SBATCH --mem 60G
#SBATCH -o /scratch/rjp5nc/busco.out # Standard output
#SBATCH -e /scratch/rjp5nc/err/busco.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

#conda create -n busco_env -c bioconda -c conda-forge busco
conda activate busco_env
# apptainer pull docker://ezlabgva/busco:v5.4.7_cv1
module load apptainer


# Move to data directory
cd /scratch/rjp5nc/HMW/HMWDNAElvis3/

# Run Busco
apptainer run /home/csm6hg/sifs/busco_v5.4.7_cv1.sif \
busco \
-i /scratch/rjp5nc/HMW/HMWDNAElvis3/hifiasm_out/assembly.hap2.p_ctg.fa \
-c 6 \
--out_path /scratch/rjp5nc/HMW/HMWDNAElvis3/ \
-l arthropoda_odb10 \
-o dap_busco \
-m genome
