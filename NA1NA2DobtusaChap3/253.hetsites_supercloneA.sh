
module load bcftools
VCF=/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/trimmed10bp_filtered_Gilmer_12.vcf.gz
OUT=/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/trimmed10bp_filtered_hetregions_Gilmer_12.subset_regions.vcf.gz

# 1) Make a 1-based regions file (CHR START END)
cat > include_regions.bed <<'EOF'
JAACYE010000003.1 2500000 2599999
JAACYE010000003.1 2600000 2699999
JAACYE010000007.1 1100000 1199999
JAACYE010000007.1 1200000 1299999
JAACYE010000008.1 1800000 1899999
JAACYE010000008.1 1900000 1999999
JAACYE010000011.1 1400000 1499999
JAACYE010000011.1 1500000 1599999
JAACYE010000012.1 1100000 1199999
EOF

# 2) Make sure the VCF is indexed (TBI or CSI). Skip if already present.
bcftools index -f $VCF
bcftools view -h "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/trimmed10bp_filtered_Gilmer_12.vcf.gz" | grep '^##contig' | cut -d= -f3 | tr -d '>' | head

# 1) Remove any Windows CRLF and blank lines
sed -i 's/\r$//' include_regions.bed
sed -i '/^[[:space:]]*$/d' include_regions.bed

# 2) Force tab-delimited 3-column BED (CHR START END)
#    (This rewrites the file with exactly three tab-separated fields per line.)
awk 'NF>=3{print $1"\t"$2"\t"$3}' include_regions.bed > include_regions.fixed.bed && mv include_regions.fixed.bed include_regions.bed

# 3) Validate formatting (should show exactly 3 columns)
awk 'NF!=3{print "Bad line " NR ": " $0}' include_regions.bed

# 4) Optional: visualize hidden characters to catch issues
# cat -A include_regions.bed | nl | head

# 5) Re-run extraction
VCF=/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/trimmed10bp_filtered_Gilmer_12.vcf.gz
OUT=/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/trimmed10bp_filtered_hetregions_Gilmer_12.subset_regions.vcf.gz

bcftools view -R include_regions.bed -Oz -o "$OUT" "$VCF"
bcftools index -t "$OUT"



module load gcc/11.4.0  openmpi/4.1.4 icu R/4.3.1
R

library(SNPRelate)
library(SeqArray)

# Input and output
vcf.fn <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/trimmed10bp_filtered_hetregions_Gilmer_12.subset_regions.vcf.gz"
gds.fn <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/trimmed10bp_filtered_Gilmer_12.subset_regions.gds"


# vcf.fn <- "/scratch/rjp5nc/UK2022_2024/allsites_mito/usdobtusa_mito_allsites.haploid.vcf.gz"
# gds.fn <- "/scratch/rjp5nc/UK2022_2024/allsites_mito/usdobtusa_mito_allsites.haploid.gds"


# Option 1: Using SeqArray (recommended for newer VCFs)
seqVCF2GDS(vcf.fn, gds.fn, storage.option="LZMA_RA", verbose=TRUE)

# Option 2: Using SNPRelate (simpler if you don’t need all annotations)
# snpgdsVCF2GDS(vcf.fn, gds.fn, method="copy.num.of.ref", ignore.chr.prefix=TRUE)

# Check content
seqSummary(gds.fn)

# (Optional) Optimize for fast access
seqOptimize(gds.fn)
