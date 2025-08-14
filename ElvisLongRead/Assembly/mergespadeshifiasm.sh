#!/usr/bin/env bash
#
#SBATCH -J mergehifiasmspades # A single job name for the array
#SBATCH --ntasks-per-node=40 # multi core
#SBATCH -N 1 # on one node
#SBATCH -t 2-00:00 # 6 days
#SBATCH --mem 100G
#SBATCH -o /scratch/rjp5nc/Canu_error/hifiasm%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/Canu_error/hifiasm%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

#conda create -n hybridmerge -c bioconda -c conda-forge ragtag minimap2 samtools nextpolish busco quast seqkit

source $(conda info --base)/etc/profile.d/conda.sh
conda activate hybridmerge

cd /scratch/rjp5nc/HMW/mergeSpadesHifiasm/

# Paths
HIFIASM=/scratch/rjp5nc/Reference_genomes/post_kraken/assembly.hap2_onlydaps.fasta
SPADES=/scratch/rjp5nc/spades/spades_output/contigs.fasta
HIFI=/scratch/rjp5nc/HMW/HMWDNAElvis3/m84128_250121_222443_s2.hifi_reads.bc2104.fq.gz
R1=/scratch/rjp5nc/HMW/shortreadElvis/merged_R1.fq.gz
R2=/scratch/rjp5nc/HMW/shortreadElvis/merged_R2.fq.gz

THREADS=$SLURM_CPUS_PER_TASK

# 1. Filter SPAdes contigs >= 1kb and merge with Hifiasm assembly using RagTag
module load seqkit ragtag minimap2 samtools nextpolish quast

seqkit seq -m 1000 "$SPADES" > spades.1kb.fasta
ragtag.py patch -o ragtag_patch "$HIFIASM" spades.1kb.fasta
MINPUT=ragtag_patch/ragtag.patch.fasta

# 2. Map HiFi reads
minimap2 -t $THREADS -x map-hifi "$MINPUT" "$HIFI" | samtools sort -@ 8 -o hifi.bam
samtools index hifi.bam

# 3. Create config for HiFi polishing
cat > np_pb.cfg <<EOF
[General]
job_type=local
task=polish
genome=$MINPUT
genome_size=150m
threads=$THREADS

[sgs_option]
sgs_bam=

[lgs_option]
lgs_bam=hifi.bam
EOF

# 4. Run NextPolish (HiFi round)
nextPolish np_pb.cfg > np_round1.fasta

# 5. Map Illumina reads to HiFi-polished assembly
minimap2 -t $THREADS -ax sr np_round1.fasta "$R1" "$R2" | samtools sort -@ 8 -o ilmn.bam
samtools index ilmn.bam

# 6. Create config for Illumina polishing
cat > np_ilmn.cfg <<EOF
[General]
job_type=local
task=polish
genome=np_round1.fasta
genome_size=150m
threads=$THREADS

[sgs_option]
sgs_bam=ilmn.bam

[lgs_option]
lgs_bam=
EOF

# 7. Run NextPolish (Illumina round)
nextPolish np_ilmn.cfg > np_round2.fasta

# 8. QC with QUAST
quast -t $THREADS -o quast_np np_round2.fasta