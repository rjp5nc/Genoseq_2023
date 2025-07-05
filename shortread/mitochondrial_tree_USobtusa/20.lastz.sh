#!/usr/bin/env bash

#SBATCH -J runFASTQC # A single job name for the array
#SBATCH --ntasks-per-node=40 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 4:00:00 ### 15 seconds
#SBATCH --mem 120G
#SBATCH -o /scratch/rjp5nc/erroroutputs/nFlo_1.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/erroroutputs/nFlo_1.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# sbatch /home/rjp5nc/Genoseq_2023/shortread/liftover/lastz.sh
# sacct -j 761713

module load nextflow/24.10.5
module load apptainer/1.3.4

### ssdf
### test container
  #nextflow run evotools/nf-LO -profile test,singularity -with-singularity ${PWD}/nflo.sif

#ref=dambigua_mito
#ref=eudobtusa_mito
#ref=kap4Dpulex_mito
ref=usdobtusa_mito

mkdir -p "/scratch/rjp5nc/Reference_genomes/mito_reference/lastz/$ref"

  nextflow run evotools/nf-LO \
  --source /scratch/rjp5nc/Reference_genomes/mito_reference/$ref.fasta \
  --target /scratch/rjp5nc/Reference_genomes/mito_reference/eudpulex_mito.fasta \
  --distance medium \
  --aligner lastz \
  --outdir /scratch/rjp5nc/Reference_genomes/mito_reference/lastz/$ref \
  -profile singularity \
  --max_memory '60.GB' \
  --max_cpus 40 \
  -with-singularity /scratch/rjp5nc/nflo/nflo.sif