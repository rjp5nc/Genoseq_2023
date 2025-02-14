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


module load samtools varscan bcftools

# sbatch --array=1-12 ~/Genoseq_2023/run_varscan.sh
# 40646429


#SLURM_ARRAY_TASK_ID=11
chr=$( sed -n ${SLURM_ARRAY_TASK_ID}p /project/berglandlab/Karen/genomefiles/ChrScaffoldList )
echo $chr

samtools mpileup \
-r ${chr} \
--fasta-ref /project/berglandlab/daphnia_ref/totalHiCwithallbestgapclosed.fa \
/scratch/rjp5nc/UK2022_2024/allshortreads/sortedbamsdedup/*.bam | \


java -jar $EBROOTVARSCAN/VarScan.v2.4.4.jar mpileup2snp \
/dev/stdin \
--min-coverage 4 \
--min-var-freq 0.001 \
--output-vcf > \
/scratch/rjp5nc/UK2022_2024/allshortreads/chr/2022seq.${chr}.vcf



# Rename sample names in the VCF to match the BAM filenames
#bcftools reheader \
#--samples <(ls /scratch/rjp5nc/UK2022_2024/allshortreads/sortedbamsdedup/*.bam | sed 's#.*/##; s/.bam//') \
#/scratch/rjp5nc/UK2022_2024/allshortreads/chr/2022seq.${chr}.vcf > \
#/scratch/rjp5nc/UK2022_2024/allshortreads/chr/2022seq.${chr}.renamed.vcf

#mv /scratch/rjp5nc/UK2022_2024/allshortreads/chr/2022seq.${chr}.renamed.vcf \
#/scratch/rjp5nc/UK2022_2024/allshortreads/chr/2022seq.${chr}.vcf