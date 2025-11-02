ijob -A berglandlab -c10 -p standard --mem=50G

module load nextflow/23.04.1
module load apptainer/1.3.4

### build container
# mkdir /scratch/rjp5nc/nflo
#
# cd /scratch/rjp5nc/nflo
# curl -O https://raw.githubusercontent.com/evotools/nf-LO/main/singularity.def
# curl -O https://raw.githubusercontent.com/evotools/nf-LO/main/environment.yml
# singularity build nflo.sif singularity.def

### test container
  nextflow run evotools/nf-LO -profile test,singularity -with-singularity ${PWD}/nflo.sif




