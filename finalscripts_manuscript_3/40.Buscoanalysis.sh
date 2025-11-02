#!/usr/bin/env bash
#
#SBATCH -J Busco # A single job name for the array
#SBATCH --ntasks-per-node=10 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-05:00 # 5 hours
#SBATCH --mem 60G
#SBATCH -o /scratch/rjp5nc/err/busco.out # Standard output
#SBATCH -e /scratch/rjp5nc/err/busco.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

#conda create -n busco_env -c bioconda -c conda-forge busco
conda activate busco_env
# apptainer pull docker://ezlabgva/busco:v5.4.7_cv1
module load apptainer

# Move to data directory

Ref=US_obtusa_onlydaps.fa
#Ref=totalHiCwithallbestgapclosed.fa
#Ref=Daphnia_ambigua_Q001_genome.fa
#Ref=assembly.hap2_onlydaps.fasta
#Ref=us_pulex_ref_kap4.fa

REFERENCE=/scratch/rjp5nc/Reference_genomes/post_kraken/$Ref

cd /scratch/rjp5nc/UK2022_2024/buscoanalysis

# Run Busco
apptainer run /home/rjp5nc/sifs/busco_v5.4.7_cv1.sif \
busco \
-i $REFERENCE \
-c 10 \
--out_path /scratch/rjp5nc/UK2022_2024/buscoanalysis \
-l arthropoda_odb10 \
-o $Ref \
-m genome




cd /scratch/rjp5nc/UK2022_2024/buscoanalysis/

# Directory with BUSCO GFF files
GFF_DIR=/scratch/rjp5nc/UK2022_2024/buscoanalysis/${Ref}/run_arthropoda_odb10/busco_sequences/single_copy_busco_sequences

# Output directory for BED files

mkdir -p /scratch/rjp5nc/UK2022_2024/buscoanalysis/beds/$Ref
BED_OUT=/scratch/rjp5nc/UK2022_2024/buscoanalysis/beds/$Ref

# Loop over GFFs and convert to BED
for GFF in $GFF_DIR/*.gff; do
    BASENAME=$(basename $GFF .gff)
    awk '$3=="gene" {print $1"\t"$4-1"\t"$5"\t"$9"\t0\t"$7}' $GFF > $BED_OUT/$BASENAME.bed
done

# Combine all BED files into one
cat $BED_OUT/*.bed > /scratch/rjp5nc/UK2022_2024/buscoanalysis/${Ref}_BUSCO_combined.bed

# Optional: sort by contig and start position
sort -k1,1 -k2,2n /scratch/rjp5nc/UK2022_2024/buscoanalysis/${Ref}_BUSCO_combined.bed > /scratch/rjp5nc/UK2022_2024/buscoanalysis/${Ref}_BUSCO_combined.sorted.bed
rm /scratch/rjp5nc/UK2022_2024/buscoanalysis/${Ref}_BUSCO_combined.bed


