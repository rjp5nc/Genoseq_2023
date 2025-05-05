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




#assembly.hap2_onlydaps.fasta


### dmel dsim
nextflow run evotools/nf-LO \
--source /scratch/rjp5nc/Reference_genomes/post_kraken/us_pulex_ref_kap4.fa \
--target /scratch/rjp5nc/dmel_dsim_TSP/refGenome/D_melanogaster_r6.12.Muller.fasta \
--distance medium \
--aligner lastz \
--outdir /scratch/rjp5nc/nFlo_lastz \
-profile singularity \
--max_memory '50.GB' \
--max_cpus 16 \
-with-singularity ${PWD}/nflo.sif
