#!/usr/bin/env bash
#SBATCH -J pf_combine
#SBATCH --ntasks-per-node=4
#SBATCH -N 1
#SBATCH -t 5:00:00
#SBATCH --mem 16G
#SBATCH -o /scratch/rjp5nc/err/pf_combine.%j.out
#SBATCH -e /scratch/rjp5nc/err/pf_combine.%j.err
#SBATCH -p standard
#SBATCH --account berglandlab_standard

set -euo pipefail
base=/scratch/rjp5nc/ssu_analysis/phyloflash
out=$base/combined
mkdir -p "$out"
cd "$base"

# ---------------------------------------------------------------
# 1. Extract each sample's archive
# ---------------------------------------------------------------
: > "$out/missing_samples.txt"
for d in */; do
  s=${d%/}
  [[ "$s" == "combined" ]] && continue
  tgz="$s/$s.phyloFlash.tar.gz"
  if [[ -s "$tgz" ]]; then
    [[ -s "$s/$s.phyloFlash.NTUabundance.csv" ]] || tar -xzf "$tgz" -C "$s"
  else
    echo "$s" >> "$out/missing_samples.txt"
  fi
done

# ---------------------------------------------------------------
# 2. Long table: sample, category, taxonomy, lowest rank, reads
#    (chloroplast and mitochondria separated from bacteria)
# ---------------------------------------------------------------
echo -e "sample\tcategory\ttaxonomy\tlowest_rank\treads" > "$out/ntu_long.tsv"
for f in */*.phyloFlash.NTUabundance.csv; do
  s=$(basename "$f" .phyloFlash.NTUabundance.csv)
  awk -v s="$s" '
    {
      n = $NF; sub(/,[^,]*$/, "", $0); tax = $0
      if (n !~ /^[0-9]+$/) next                     # skip header
      if (tax ~ /Chloroplast/)        cat = "Chloroplast"
      else if (tax ~ /Mitochondria/)  cat = "Mitochondria"
      else if (tax ~ /^Bacteria/)     cat = "Bacteria"
      else if (tax ~ /^Archaea/)      cat = "Archaea"
      else if (tax ~ /^Eukaryota/)    cat = "Eukaryota"
      else                            cat = "Unassigned"
      k = split(tax, p, ";"); low = p[k]
      # label unresolved genera with their parent lineage
      if (low ~ /^\(/ && k > 1) { lab = p[1]; for (i=2; i<=k; i++) lab = lab ";" p[i]; low = lab }
      print s "\t" cat "\t" tax "\t" low "\t" n
    }' FS="," "$f" >> "$out/ntu_long.tsv"
done

# ---------------------------------------------------------------
# 3. Per-sample summary
# ---------------------------------------------------------------
awk -F'\t' '
  NR > 1 { t[$1] += $5; c[$1, $2] += $5; s[$1] = 1 }
  END {
    print "sample\ttotal_SSU_reads\teukaryota\tbacteria\tarchaea\tchloroplast\tmitochondria\tunassigned\tbacteria_fraction"
    for (x in s) {
      b = c[x, "Bacteria"] + 0
      printf "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.5f\n", x, t[x],
        c[x,"Eukaryota"]+0, b, c[x,"Archaea"]+0, c[x,"Chloroplast"]+0,
        c[x,"Mitochondria"]+0, c[x,"Unassigned"]+0, (t[x] ? b / t[x] : 0)
    }
  }' "$out/ntu_long.tsv" | (read -r h; echo "$h"; sort) > "$out/sample_summary.tsv"

# ---------------------------------------------------------------
# 4. Bacteria genus x sample: counts and relative abundance
# ---------------------------------------------------------------
awk -F'\t' -v outc="$out/bacteria_genus_counts.tsv" -v outr="$out/bacteria_genus_relabund.tsv" '
  NR > 1 && $2 == "Bacteria" {
    n[$4, $1] += $5; g[$4] += $5; bt[$1] += $5
  }
  NR > 1 { samp[$1] = 1 }
  END {
    ns = 0; for (x in samp) S[++ns] = x
    # sort samples alphabetically
    for (i = 1; i <= ns; i++) for (j = i + 1; j <= ns; j++) if (S[j] < S[i]) { t = S[i]; S[i] = S[j]; S[j] = t }
    hdr = "genus"; for (i = 1; i <= ns; i++) hdr = hdr "\t" S[i]
    print hdr > outc; print hdr > outr
    cmd = "sort -t\"\t\" -k2,2nr"
    for (x in g) print x "\t" g[x] | cmd
    close(cmd)
  }' "$out/ntu_long.tsv" > "$out/.genus_order.tmp"

awk -F'\t' -v outc="$out/bacteria_genus_counts.tsv" -v outr="$out/bacteria_genus_relabund.tsv" '
  FNR == NR { order[++no] = $1; next }                       # genus order (by total reads)
  FNR == 1  { next }
  $2 == "Bacteria" { n[$4, $1] += $5; bt[$1] += $5 }
  { samp[$1] = 1 }
  END {
    ns = 0; for (x in samp) S[++ns] = x
    for (i = 1; i <= ns; i++) for (j = i + 1; j <= ns; j++) if (S[j] < S[i]) { t = S[i]; S[i] = S[j]; S[j] = t }
    for (k = 1; k <= no; k++) {
      gn = order[k]; lc = gn; lr = gn
      for (i = 1; i <= ns; i++) {
        v = n[gn, S[i]] + 0
        lc = lc "\t" v
        lr = lr "\t" (bt[S[i]] ? sprintf("%.6f", v / bt[S[i]]) : 0)
      }
      print lc >> outc; print lr >> outr
    }
  }' "$out/.genus_order.tmp" "$out/ntu_long.tsv"
rm -f "$out/.genus_order.tmp"

# ---------------------------------------------------------------
# 5. All assembled full-length SSU genes + host 18S check
# ---------------------------------------------------------------
echo -e "sample\tOTU\tread_cov\tcoverage\tdbHit\ttaxonomy\tpid\talnlen\tevalue" > "$out/assembled_SSU_all.tsv"
for f in */*.phyloFlash.extractedSSUclassifications.csv; do
  s=$(basename "$f" .phyloFlash.extractedSSUclassifications.csv)
  awk -F, -v s="$s" 'NR > 1 && NF >= 5 {
      line = s; for (i = 1; i <= NF; i++) line = line "\t" $i; print line }' "$f" \
    >> "$out/assembled_SSU_all.tsv"
done

# Highest-coverage eukaryotic SSU per sample
{ head -1 "$out/assembled_SSU_all.tsv"
  awk -F'\t' 'NR > 1 && $6 ~ /^Eukaryota/' "$out/assembled_SSU_all.tsv" \
    | sort -t$'\t' -k1,1 -k4,4gr | awk -F'\t' '!seen[$1]++'
} > "$out/host_18S_check.tsv"

# ---------------------------------------------------------------
# 6. Run statistics (phyloFlash report) per sample
# ---------------------------------------------------------------
echo -e "sample\tmetric\tvalue" > "$out/run_stats_long.tsv"
for f in */*.phyloFlash.report.csv; do
  s=$(basename "$f" .phyloFlash.report.csv)
  awk -F, -v s="$s" 'NF >= 2 && $1 != "" { print s "\t" $1 "\t" $2 }' "$f" >> "$out/run_stats_long.tsv"
done

# ---------------------------------------------------------------
echo "Samples combined:  $(($(wc -l < "$out/sample_summary.tsv") - 1))"
echo "Missing/failed:    $(wc -l < "$out/missing_samples.txt")"
echo "Bacterial genera:  $(($(wc -l < "$out/bacteria_genus_counts.tsv") - 1))"
echo "Outputs in:        $out"