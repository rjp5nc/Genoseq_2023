
NAobtusaMito.gff


input_gff="NAobtusaMito.gff"
output_bed="output.bed"

awk '
BEGIN { FS=OFS="\t" }
# Skip comments and empty lines
/^#/ || NF==0 { next }
{
  # Skip region lines
  if ($3 == "region") next;

  chrom = $1;
  start = $4 - 1;    # GFF 1-based to BED 0-based start
  end = $5;
  score = $6 == "." ? "0" : $6;
  strand = ($7 == "+" || $7 == "-") ? $7 : ".";

  # Parse attributes field ($9) for Name or ID
  name = "";
  n = split($9, attrs, ";");
  for (i=1; i<=n; i++) {
    if (attrs[i] ~ /^Name=/) {
      sub(/^Name=/, "", attrs[i]);
      name = attrs[i];
      break;
    }
  }
  if (name == "") {
    for (i=1; i<=n; i++) {
      if (attrs[i] ~ /^ID=/) {
        sub(/^ID=/, "", attrs[i]);
        name = attrs[i];
        break;
      }
    }
  }
  if (name == "") {
    name = $3;
  }

  print chrom, start, end, name, score, strand;
}
' "$input_gff" > "$output_bed"

echo "BED file created: $output_bed"