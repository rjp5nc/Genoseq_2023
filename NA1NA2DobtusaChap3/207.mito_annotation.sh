#ijob -A berglandlab -c10 -p standard --mem=40G

module load bcftools

cd /scratch/rjp5nc/Reference_genomes/mito_reference/


GFF=/scratch/rjp5nc/Reference_genomes/mito_reference/usdobtusa_annotation/usdobtusa.gff
INVCF=/scratch/rjp5nc/UK2022_2024/allsites_mito/usdobtusa_mito_allsites.haploid.vcf.gz
TSV=/scratch/rjp5nc/Reference_genomes/mito_reference/usdobtusa_annotation/usdobtusa.annot.tsv
TSV_GZ=${TSV}.gz
OUTVCF=/scratch/rjp5nc/UK2022_2024/allsites_mito/usdobtusa_mito_allsites.annotated.vcf.gz

# 1) Convert GFF → clean TSV: CHROM FROM TO source type score strand phase ID Name gene_id Parent
awk -F'\t' 'BEGIN{OFS="\t"}
  $0 !~ /^#/ {
    # parse attributes (col 9) into key->value
    split($9, kvs, ";"); id=""; name=""; gid=""; par=""
    for(i in kvs){
      split(kvs[i], kv, "=")
      key=kv[1]; val=kv[2]
      if(key=="ID") id=val
      else if(key=="Name") name=val
      else if(key=="gene_id") gid=val
      else if(key=="Parent") par=val
    }
    # print: 1=CHROM, 4=FROM, 5=TO, 2=source, 3=type, 6=score, 7=strand, 8=phase, attrs...
    print $1,$4,$5,$2,$3,$6,$7,$8,id,name,gid,par
  }' "$GFF" > "$TSV"

# 2) Sort, bgzip, index as a generic tabix file
sort -k1,1 -k2,2n "$TSV" | bgzip -c > "$TSV_GZ"
tabix -s 1 -b 2 -e 3 "$TSV_GZ"

# 3) INFO header for the fields we’ll add
cat > /scratch/rjp5nc/Reference_genomes/mito_reference/usdobtusa_annotation/gff_info_hdr.txt <<'EOF'
##INFO=<ID=gff_source,Number=1,Type=String,Description="GFF source (col2)">
##INFO=<ID=gff_type,Number=1,Type=String,Description="GFF feature type (col3)">
##INFO=<ID=gff_score,Number=1,Type=String,Description="GFF score (col6)">
##INFO=<ID=gff_strand,Number=1,Type=String,Description="GFF strand (col7)">
##INFO=<ID=gff_phase,Number=1,Type=String,Description="GFF phase (col8)">
##INFO=<ID=ID,Number=1,Type=String,Description="GFF attribute: ID">
##INFO=<ID=Name,Number=1,Type=String,Description="GFF attribute: Name">
##INFO=<ID=gene_id,Number=1,Type=String,Description="GFF attribute: gene_id">
##INFO=<ID=Parent,Number=.,Type=String,Description="GFF attribute: Parent (may be multiple)">
EOF

# 4) Annotate from the TSV (columns are 1-based in -c)
bcftools annotate \
  -a "$TSV_GZ" \
  -c CHROM,FROM,TO,INFO/gff_source:=4,INFO/gff_type:=5,INFO/gff_score:=6,INFO/gff_strand:=7,INFO/gff_phase:=8,INFO/ID:=9,INFO/Name:=10,INFO/gene_id:=11,INFO/Parent:=12 \
  -h /scratch/rjp5nc/Reference_genomes/mito_reference/usdobtusa_annotation/gff_info_hdr.txt \
  -m + \
  -Oz -o "$OUTVCF" \
  "$INVCF"

# 5) Index the new VCF
tabix -p vcf "$OUTVCF"



module load htslib


RAWVCF=/scratch/rjp5nc/UK2022_2024/allsites_mito/usdobtusa_mito_allsites.haploid.vcf.gz
TSV_GZ=/scratch/rjp5nc/Reference_genomes/mito_reference/usdobtusa_annotation/usdobtusa.annot.tsv.gz
OUTVCF=/scratch/rjp5nc/UK2022_2024/allsites_mito/usdobtusa_mito_allsites.annot2.vcf.gz

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



BADVCF=/scratch/rjp5nc/UK2022_2024/allsites_mito/usdobtusa_mito_allsites.annot2.vcf.gz
CLEANVCF=/scratch/rjp5nc/UK2022_2024/allsites_mito/usdobtusa_mito_allsites.annot2.clean.vcf.gz

# extract header, drop the empty-ID INFO line
bcftools view -h "$BADVCF" | grep -v '^##INFO=<ID=,' > /scratch/rjp5nc/tmp.clean.header.vcf

# rebuild: cleaned header + records
{ cat /scratch/rjp5nc/tmp.clean.header.vcf; bcftools view -H "$BADVCF"; } \
  | bgzip -c > "$CLEANVCF"

# index
bcftools index -t "$CLEANVCF"

# sanity check: no blank INFO IDs should remain
bcftools view -h "$CLEANVCF" | grep -E '^##INFO=<ID='
