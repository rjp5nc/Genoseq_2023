#!/usr/bin/env bash

#SBATCH -J rename # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 6-0:00:00 ### 15 seconds
#SBATCH --mem 60G
#SBATCH -o /scratch/rjp5nc/erroroutputs/beast.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/erroroutputs/beast.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab
#SBATCH --mail-type=END               # Send email at job completion
#SBATCH --mail-user=rjp5nc@virginia.edu    # Email address for notifications

module load bcftools

#bcftools annotate --rename-chrs rename_chroms.txt input.vcf.gz -Oz -o renamed.vcf.gz

cd /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/

cp trimmed10bp_allsites_usobtusa.vcf.gz trimmed10bp_allsites_usobtusa2.vcf.gz

bcftools annotate --rename-chrs contigmap.txt trimmed10bp_allsites_usobtusa.vcf.gz -Oz -o trimmed10bp_allsites_usobtusa_renamed.vcf.gz
bcftools index trimmed10bp_allsites_usobtusa_renamed.vcf.gz

# awk -v MAP="contigmap.txt" '
#   BEGIN{
#     # read map; split on ANY whitespace
#     while ((getline line < MAP) > 0) {
#       if (line ~ /^[[:space:]]*$/) continue
#       n = split(line, a, /[[:space:]]+/)
#       if (n >= 2) m[a[1]] = a[2]
#     }
#     close(MAP)
#   }
#   # For FASTA
#   /^>/{
#     hdr = substr($0,2)                 # drop ">"
#     id  = hdr
#     if (match(hdr, /^[^[:space:]]+/))  # first token = true ID
#       id = substr(hdr, RSTART, RLENGTH)
#     rest = substr(hdr, length(id)+1)   # preserve description
#     if (id in m) id = m[id]
#     print ">" id rest
#     next
#   }
#   { print }                            # sequence lines
# ' /scratch/rjp5nc/Reference_genomes/post_kraken/US_obtusa_onlydaps.fa > US_obtusa_onlydaps.renamed.fa
# cat US_obtusa_onlydaps.renamed.fa | head -n 5


# awk '/^>/{sub(/_LG[0-9]+/,""); print; next} {print}' US_obtusa_onlydaps.renamed.fa > US_obtusa_onlydaps.nolgs.fa
# samtools faidx US_obtusa_onlydaps.renamed.fa


# conda create -n gfftools -c bioconda -c conda-forge gffread samtools bcftools bedtools
# conda activate gfftools
# gffread Daphnia_obtusa_FS6_genome.gtf -E -o Daphnia_obtusa_FS6_genome.gff3

# awk 'BEGIN{FS=OFS="\t"}
# /^#/ {print; next}
# {sub(/_LG[0-9]+$/, "", $1); print}' \
# Daphnia_obtusa_FS6_genome.gff3 > Daphnia_obtusa_FS6_genome.nolgs.gff3

# ( grep '^#' Daphnia_obtusa_FS6_genome.nolgs.gff3
#   grep -v '^#' Daphnia_obtusa_FS6_genome.nolgs.gff3 \
#     | sort -t $'\t' -k1,1 -k4,4n -k5,5n
# ) | bgzip > Daphnia_obtusa_FS6_genome.nolgs.sorted.gff3.gz

# Index it
# tabix -p gff Daphnia_obtusa_FS6_genome.nolgs.sorted.gff3.gz

bcftools csq \
  -f US_obtusa_onlydaps.renamed.fa \
  -g Daphnia_obtusa_FS6_genome.nolgs.sorted.gff3.gz \
  -Oz -o trimmed10bp_allsites_usobtusa_csq.vcf.gz \
  trimmed10bp_allsites_usobtusa_renamed.vcf.gz

bcftools index trimmed10bp_allsites_usobtusa_csq.vcf.gz