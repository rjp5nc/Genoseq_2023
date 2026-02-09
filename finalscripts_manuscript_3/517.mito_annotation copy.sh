#ijob -A berglandlab -c10 -p standard --mem=40G

module load bcftools

cd /scratch/rjp5nc/Reference_genomes/mito_reference/



module load htslib


RAWVCF=/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/usdobtusa_mito_allsites_all.haploid.vcf.gz
TSV_GZ=/scratch/rjp5nc/Reference_genomes/mito_reference/usdobtusa_annotation/usdobtusa.annot.tsv.gz
OUTVCF=/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/usdobtusa_mito_allsites_all.haploid.annotated.vcf.gz

# New header with SAFE names for attributes (prefixed with gff_)
cat > /scratch/rjp5nc/Reference_genomes/mito_reference/usdobtusa_annotation/gff_info_hdr.safe.txt <<'EOF'
##INFO=<ID=gff_source,Number=1,Type=String,Description="GFF source (col2)">
##INFO=<ID=gff_type,Number=1,Type=String,Description="GFF feature type (col3)">
##INFO=<ID=gff_score,Number=1,Type=String,Description="GFF score (col6)">
##INFO=<ID=gff_strand,Number=1,Type=String,Description="GFF strand (col7)">
##INFO=<ID=gff_phase,Number=1,Type=String,Description="GFF phase (col8)">
##INFO=<ID=gff_ID,Number=1,Type=String,Description="GFF attribute: ID">
##INFO=<ID=gff_Name,Number=1,Type=String,Description="GFF attribute: Name">
##INFO=<ID=gff_gene_id,Number=1,Type=String,Description="GFF attribute: gene_id">
##INFO=<ID=gff_Parent,Number=.,Type=String,Description="GFF attribute: Parent (may be multiple)">
EOF

# Re-annotate starting from the UNANNOTATED VCF
bcftools annotate \
  -a "$TSV_GZ" \
  -c CHROM,FROM,TO,INFO/gff_source:=4,INFO/gff_type:=5,INFO/gff_score:=6,INFO/gff_strand:=7,INFO/gff_phase:=8,INFO/gff_ID:=9,INFO/gff_Name:=10,INFO/gff_gene_id:=11,INFO/gff_Parent:=12 \
  -h /scratch/rjp5nc/Reference_genomes/mito_reference/usdobtusa_annotation/gff_info_hdr.safe.txt \
  -m + \
  -Oz -o "$OUTVCF" \
  "$RAWVCF"

tabix -p vcf "$OUTVCF"

# Quick sanity check: ensure no blank INFO IDs exist
bcftools view -h "$OUTVCF" | grep -E '^##INFO=<ID='
# You should NOT see a line like: ##INFO=<ID=, ...>



BADVCF=/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/usdobtusa_mito_allsites_all.haploid.annotated.vcf.gz
CLEANVCF=/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/usdobtusa_mito_allsites_all.haploid.annotated2.vcf.gz

# extract header, drop the empty-ID INFO line
bcftools view -h "$BADVCF" | grep -v '^##INFO=<ID=,' > /scratch/rjp5nc/tmp.clean.header.vcf

# rebuild: cleaned header + records
{ cat /scratch/rjp5nc/tmp.clean.header.vcf; bcftools view -H "$BADVCF"; } \
  | bgzip -c > "$CLEANVCF"

# index
bcftools index -t "$CLEANVCF"

# sanity check: no blank INFO IDs should remain
bcftools view -h "$CLEANVCF" | grep -E '^##INFO=<ID='
