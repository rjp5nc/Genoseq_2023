#!/usr/bin/env bash
#
#SBATCH -J gatk_chrom # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-10:00:00 # 8 hours
#SBATCH --mem 25G
#SBATCH -o /scratch/rjp5nc/err/gatk.chrom.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/err/gatk.chrom.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab
#SBATCH --array=1-521%60
#SBATCH --mail-type=END               # Send email at job completion
#SBATCH --mail-user=rjp5nc@virginia.edu    # Email address for notifications

# grep '^>' your_file.fasta | sed 's/^>//'
# grep '^>' dambigua_mito.fasta | sed 's/^>//' > /scratch/rjp5nc/UK2022_2024/dambigua_mito_chr.txt
# grep '^>' eudpulex_mito.fasta | sed 's/^>//' > /scratch/rjp5nc/UK2022_2024/eudpulex_mito_chr.txt
# grep '^>' eudobtusa_mito.fasta | sed 's/^>//' > /scratch/rjp5nc/UK2022_2024/eudobtusa_mito_chr.txt
# grep '^>' kap4Dpulex_mito.fasta | sed 's/^>//' > /scratch/rjp5nc/UK2022_2024/kap4Dpulex_mito_chr.txt
# grep '^>' usdobtusa_mito.fasta | sed 's/^>//' > /scratch/rjp5nc/UK2022_2024/usdobtusa_mito_chr.txt

#gatk CreateSequenceDictionary -R dambigua_mito.fasta
#gatk CreateSequenceDictionary -R eudpulex_mito.fasta
#gatk CreateSequenceDictionary -R eudobtusa_mito.fasta
#gatk CreateSequenceDictionary -R kap4Dpulex_mito.fasta
#gatk CreateSequenceDictionary -R usdobtusa_mito.fasta

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

#finished param10

parameterFile=/scratch/rjp5nc/UK2022_2024/touseforDBI_mito_fullref.csv
wd="/scratch/rjp5nc/UK2022_2024/mitogvcf"

#dos2unix "$parameterFile"
  
#awk -F',' 'NR>1 {print $1 "\t/scratch/rjp5nc/UK2022_2024/mitogvcf/gvcf/" $6 "/" $1 ".g.vcf.gz"}' /scratch/rjp5nc/UK2022_2024/touseforDBI_mito_fullref.csv > /scratch/rjp5nc/UK2022_2024/sample_map.txt


# Extract sample name

#changed these for pulex file
samp=$(awk -F',' -v task_id="$SLURM_ARRAY_TASK_ID" 'NR == task_id {print $1}' "$parameterFile")
ref=$(awk -F',' -v task_id="$SLURM_ARRAY_TASK_ID" 'NR == task_id {print $5}' "$parameterFile")
folder=$(awk -F',' -v task_id="$SLURM_ARRAY_TASK_ID" 'NR == task_id {print $6}' "$parameterFile")

samp=$(echo $samp | tr -d '\r')
ref=$(echo $ref | tr -d '\r')
folder=$(echo $folder | tr -d '\r')

echo "Haplotype calling -" "Sample:" $SLURM_ARRAY_TASK_ID
echo "Chromosome:" ${folder}

# Create folder for chromosome
if [[ -d "${wd}/gvcf/${folder}" ]]
then
	echo "Working chromosome folder exist"
	echo "lets move on"
	date
else
	echo "folder doesnt exist. lets fix that"
	mkdir ${wd}/gvcf/${folder}
	date
fi


#gatk CreateSequenceDictionary -R "$ref"
#samtools faidx "$ref"


# Haplotype Calling
#gatk HaplotypeCaller \
#-R $ref \
#-I /scratch/rjp5nc/UK2022_2024/final_mitobam_rg2/${samp}finalmap_RG.bam \
#-O ${wd}/gvcf/${folder}/${samp}.g.vcf \
#-L /scratch/rjp5nc/Reference_genomes/mito_reference/${folder}.bed \
#-ERC GVCF

# Bgzip
module load gcc/14.2.0 htslib/1.17

# Compress and index with Tabix
bgzip ${wd}/gvcf/${folder}/${samp}.g.vcf
tabix -p vcf ${wd}/gvcf/${folder}/${samp}.g.vcf.gz

# Finish
echo "Complete -" "Sample:" $SLURM_ARRAY_TASK_ID