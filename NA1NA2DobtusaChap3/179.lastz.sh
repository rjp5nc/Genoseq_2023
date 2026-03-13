#!/usr/bin/env bash

#SBATCH -J lastz # A single job name for the array
#SBATCH --ntasks-per-node=40 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 6-00:00:00 ### 15 seconds
#SBATCH --mem 120G
#SBATCH -o /scratch/rjp5nc/erroroutputs/nFlo_1.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/erroroutputs/nFlo_1.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# sbatch /home/rjp5nc/Genoseq_2023/shortread/liftover/lastz.sh
# sacct -j 761713

module load nextflow/25.04.6
module load apptainer/1.3.4

### ssdf
### test container
  #nextflow run evotools/nf-LO -profile test,singularity -with-singularity ${PWD}/nflo.sif

#mkdir -p /scratch/rjp5nc/lastz/eu_obtusa


  nextflow run evotools/nf-LO \
  --source /scratch/rjp5nc/Reference_genomes/post_kraken/assembly.hap2_onlydaps.fasta \
  --target /scratch/rjp5nc/Reference_genomes/post_kraken/US_obtusa_onlydaps.fa \
  --distance medium \
  --aligner lastz \
  --outdir /scratch/rjp5nc/lastz/eu_obtusa_to_usobtusa \
  -profile singularity \
  --max_memory '60.GB' \
  --max_cpus 40 \
  -with-singularity /scratch/rjp5nc/nflo/nflo.sif