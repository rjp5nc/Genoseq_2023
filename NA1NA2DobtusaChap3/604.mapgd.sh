#!/usr/bin/env bash
#
#SBATCH -J pixy # A single job name for the array
#SBATCH --ntasks-per-node=20 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-92:00  ### 48 hours
#SBATCH --mem 50G
#SBATCH -o /scratch/rjp5nc/err/pixy.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/err/pixy.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab


VCF=/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/trimmed10bp_allsites_usobtusa_renamed_annotated.vcf.gz
OUTDIR=/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/mapgd_relatedness
mkdir -p "$OUTDIR"


bcftools view -h "$VCF" | grep -E '^##FORMAT=<ID=(GL|PL)\b' || echo "No GL/PL header found"

# Look at FORMAT + first few samples for a handful of variant lines
bcftools view -H "$VCF" | head -n 5 | cut -f1-10
bcftools query -r "$(bcftools query -f '%CHROM\t%POS\n' "$VCF" | head -n1 | awk '{print $1":"$2"-"$2}')" \
  -f '%CHROM\t%POS\t%REF\t%ALT\t%FORMAT[\t%GL][\t%PL]\n' "$VCF" | head






MAPGD=/scratch/rjp5nc/mapgd/MAPGD/bin/mapgd
VCF=/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/trimmed10bp_allsites_usobtusa_renamed_annotated.vcf.gz
OUTDIR=/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/mapgd
PREFIX=usobtusa

mkdir -p "$OUTDIR"



cd $OUTDIR

/scratch/rjp5nc/mapgd/MAPGD/bin/mapgd read --help
/scratch/rjp5nc/mapgd/MAPGD/bin/mapgd allele $VCF 
