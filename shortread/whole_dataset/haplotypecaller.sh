#!/usr/bin/env bash
#
#SBATCH -J gatk_chrom # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-12:00:00 # 8 hours
#SBATCH --mem 25G
#SBATCH -o /scratch/rjp5nc/err/gatk.chrom.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/err/gatk.chrom.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab
#SBATCH --array=1-672%50


# grep '^>' your_file.fasta | sed 's/^>//'
# grep '^>' assembly.hap2_onlydaps.fasta | sed 's/^>//' > /scratch/rjp5nc/UK2022_2024/assembly.hap2_onlydaps_chr.txt
# grep '^>' Daphnia_ambigua_Q001_genome.fa | sed 's/^>//' > /scratch/rjp5nc/UK2022_2024/ambigua_chr.txt
# grep '^>' eu_pulex_totalHiCwithallbestgapclosed.clean.fa | sed 's/^>//' > /scratch/rjp5nc/UK2022_2024/eu_pulex_chr.txt
# grep '^>' US_obtusa_onlydaps.fa | sed 's/^>//' > /scratch/rjp5nc/UK2022_2024/us_obtusa_chr.txt
# grep '^>' us_pulex_ref_kap4.fa | sed 's/^>//' > /scratch/rjp5nc/UK2022_2024/us_pulex_chr.txt
# grep '^>' eu_pulex_totalHiCwithallbestgapclosed.fa | sed 's/^>//' > /scratch/rjp5nc/UK2022_2024/eu_pulex_all_chr.txt

#gatk CreateSequenceDictionary -R Daphnia_ambigua_Q001_genome.fa
#gatk CreateSequenceDictionary -R totalHiCwithallbestgapclosed.fa

#gatk CreateSequenceDictionary -R assembly.hap2_onlydaps.fasta
#gatk CreateSequenceDictionary -R us_pulex_ref_kap4.fa

#cut -f1,2 eu_pulex_totalHiCwithallbestgapclosed.clean.fai | awk '{print $1"\t0\t"$2}' > your.bed

#while read -r line; do
#    chrom=$(echo "$line" | cut -f1)
#    length=$(echo "$line" | cut -f2)
#    echo -e "${chrom}\t0\t${length}" > "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/bed/${chrom}.bed"
#done < "eu_pulex_totalHiCwithallbestgapclosed.clean.fa.fai"

#grep "^>h2tg000002l" assembly.hap2_onlydaps.fasta

# Load modules
module load gatk/4.6.0.0
module load tabix/0.2.6

# Parameters     find /project/berglandlab/connor -type f -name "SRA_paramList_1"

#parameterFile=/scratch/rjp5nc/UK2022_2024/robert_paramfile.txt
# wc -l /scratch/rjp5nc/UK2022_2024/robert_paramfile.txt
# total lines : 101136

#sed -n '1,9999p' /scratch/rjp5nc/UK2022_2024/robert_paramfile.txt > /scratch/rjp5nc/UK2022_2024/param1_9999.txt
#sed -n '10000,19998p' /scratch/rjp5nc/UK2022_2024/robert_paramfile.txt > /scratch/rjp5nc/UK2022_2024/param10000_19998.txt
#sed -n '19999,29997p' /scratch/rjp5nc/UK2022_2024/robert_paramfile.txt > /scratch/rjp5nc/UK2022_2024/param3.txt
#sed -n '29998,39996p' /scratch/rjp5nc/UK2022_2024/robert_paramfile.txt > /scratch/rjp5nc/UK2022_2024/param4.txt
#sed -n '39997,49995p' /scratch/rjp5nc/UK2022_2024/robert_paramfile.txt > /scratch/rjp5nc/UK2022_2024/param5.txt
#sed -n '49996,59994p' /scratch/rjp5nc/UK2022_2024/robert_paramfile.txt > /scratch/rjp5nc/UK2022_2024/param6.txt
#sed -n '59995,69993p' /scratch/rjp5nc/UK2022_2024/robert_paramfile.txt > /scratch/rjp5nc/UK2022_2024/param7.txt
#sed -n '69994,79992p' /scratch/rjp5nc/UK2022_2024/robert_paramfile.txt > /scratch/rjp5nc/UK2022_2024/param8.txt
#sed -n '79993,89991p' /scratch/rjp5nc/UK2022_2024/robert_paramfile.txt > /scratch/rjp5nc/UK2022_2024/param9.txt
#sed -n '89992,99990p' /scratch/rjp5nc/UK2022_2024/robert_paramfile.txt > /scratch/rjp5nc/UK2022_2024/param10.txt
#sed -n '99991,101136p' /scratch/rjp5nc/UK2022_2024/robert_paramfile.txt > /scratch/rjp5nc/UK2022_2024/param11.txt

parameterFile=/scratch/rjp5nc/UK2022_2024/param1_9999.txt
wd="/scratch/rjp5nc/UK2022_2024/daphnia_phylo"

#dos2unix "$parameterFile"
  
# Extract sample name

#changed these for pulex file
id=$(awk -F',' -v task_id="$SLURM_ARRAY_TASK_ID" 'NR == task_id {print $7}' "$parameterFile")
samp=$(awk -F',' -v task_id="$SLURM_ARRAY_TASK_ID" 'NR == task_id {print $2}' "$parameterFile")
chrom=$(awk -F',' -v task_id="$SLURM_ARRAY_TASK_ID" 'NR == task_id {print $6}' "$parameterFile")
ref=$(awk -F',' -v task_id="$SLURM_ARRAY_TASK_ID" 'NR == task_id {print $5}' "$parameterFile")


wd=$(echo $wd | tr -d '\r')
chrom=$(echo $chrom | tr -d '\r')
samp=$(echo $samp | tr -d '\r')
id=$(echo $id | tr -d '\r')




echo "Haplotype calling -" "Sample:" $SLURM_ARRAY_TASK_ID
echo "Chromosome:" ${chrom}

# Create folder for chromosome
if [[ -d "${wd}/gvcf/${chrom}" ]]
then
	echo "Working chromosome folder exist"
	echo "lets move on"
	date
else
	echo "folder doesnt exist. lets fix that"
	mkdir ${wd}/gvcf/${chrom}
	date
fi


#gatk CreateSequenceDictionary -R "$ref"
#samtools faidx "$ref"


# Haplotype Calling
gatk HaplotypeCaller \
-R $ref \
-I /scratch/rjp5nc/UK2022_2024/final_bam_rg2/${samp}finalmap_RG.bam \
-O ${wd}/gvcf3/${chrom}/${samp}.${chrom}.${id}.g.vcf \
-L ${wd}/bed/${chrom}.bed \
-ERC GVCF

# Bgzip
module load gcc/14.2.0 htslib/1.17

# Compress and index with Tabix
bgzip ${wd}/gvcf3/${chrom}/${samp}.${chrom}.${id}.g.vcf
tabix -p vcf ${wd}/gvcf3/${chrom}/${samp}.${chrom}.${id}.g.vcf.gz

# Finish
echo "Complete -" "Sample:" $SLURM_ARRAY_TASK_ID