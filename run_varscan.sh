#!/usr/bin/env bash
#
#SBATCH -J download # A single job name for the array
#SBATCH --ntasks-per-node=10 # one cores
#SBATCH -N 1 # on one node
#SBATCH -t 36:00:00 ### 20 hours
#SBATCH --mem 90G
#SBATCH -o /scratch/rjp5nc/download_%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/download_%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab


module load samtools varscan

# sbatch --array=1-12 ~/Genoseq2023/variantCalling/run_varscan.sh
# 40646429


#SLURM_ARRAY_TASK_ID=11
chr=$( sed -n ${SLURM_ARRAY_TASK_ID}p /project/berglandlab/Karen/genomefiles/ChrScaffoldList )
echo $chr

samtools mpileup \
-r ${chr} \
--fasta-ref /project/berglandlab/daphnia_ref/totalHiCwithallbestgapclosed.fa \
$(ls /project/berglandlab/Robert/UKSequencing2022_2024/usftp21.novogene.com/01.RawData/Bams/sort.dedup.bam/*.bam) | bcftools call -mv -Ov -o output_variants.vcf



java -jar $EBROOTVARSCAN/VarScan.v2.4.4.jar mpileup2snp \
/dev/stdin \
--min-coverage 4 \
--min-var-freq 0.001 \
--output-vcf > \
/project/berglandlab/Robert/UKSequencing2022_2024/usftp21.novogene.com/01.RawData/Bams/2022seq.${chr}.vcf