

conda create -n sra-meta \
  python=3.10 \
  pandas \
  lxml \
  curl \
  -c conda-forge -y



IDS="/scratch/rjp5nc/rawdata/sra_ids.txt"
OUTDIR="/scratch/rjp5nc/rawdata/sra_metadata_out"
mkdir -p "$OUTDIR"
cd "$OUTDIR"

RUNINFO_CSV="runinfo.csv"
RUNINFO_TSV="runinfo.tsv"
BIOSAMPLES_TXT="biosamples.txt"
BIOSAMPLE_ATTR_LONG="biosample_attributes_long.tsv"
BIOSAMPLE_ATTR_WIDE="biosample_attributes_wide.tsv"
MERGED="sra_merged.tsv"


grep -v '^\s*$' "$IDS" \
  | sed 's/#.*$//' \
  | awk 'NF{print $1}' \
  | sort -u \
  > ids.clean.txt

wc -l ids.clean.txt
head ids.clean.txt


: > "$RUNINFO_CSV"
first=1

while read -r acc; do
  echo "[RunInfo] $acc"
  tmp="${acc}.runinfo.csv"

  if curl -fsSL \
    "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?acc=${acc}" \
    -o "$tmp"; then

    if [[ $first -eq 1 ]]; then
      cat "$tmp" >> "$RUNINFO_CSV"
      first=0
    else
      tail -n +2 "$tmp" >> "$RUNINFO_CSV"
    fi
  else
    echo "  failed: $acc"
  fi
done < ids.clean.txt


head "$RUNINFO_CSV"


python3 - <<'PY'
import csv
with open("runinfo.csv", newline='') as f, open("runinfo.tsv", "w", newline='') as g:
    r = csv.reader(f)
    w = csv.writer(g, delimiter="\t")
    for row in r:
        w.writerow(row)
print("wrote runinfo.tsv")
PY

python3 - <<'PY'
import pandas as pd

df = pd.read_csv("runinfo.tsv", sep="\t", dtype=str).fillna("")
assert "BioSample" in df.columns, "BioSample column missing"

bs = sorted(x for x in df["BioSample"] if x and x != "-")

with open("biosamples.txt","w") as w:
    for x in bs:
        w.write(x+"\n")

print("BioSamples:", len(bs))
PY

: > "$BIOSAMPLE_ATTR_LONG"


python3 - <<'PY'
import subprocess, xml.etree.ElementTree as ET

ids = [l.strip() for l in open("biosamples.txt") if l.strip()]
out = open("biosample_attributes_long.tsv", "a")

def fetch(batch):
    url = (
      "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
      "?db=biosample&id=" + ",".join(batch) + "&retmode=xml"
    )
    return subprocess.check_output(["curl","-fsSL",url], text=True)

chunk = 100
for i in range(0, len(ids), chunk):
    batch = ids[i:i+chunk]
    try:
        xmltxt = fetch(batch)
    except subprocess.CalledProcessError:
        continue
    root = ET.fromstring(xmltxt)
    for bs in root.findall(".//BioSample"):
        acc = bs.get("accession","")
        for a in bs.findall(".//Attributes/Attribute"):
            name = (a.get("attribute_name") or "").strip()
            val  = (a.text or "").strip()
            if acc and name and val:
                out.write(f"{acc}\t{name}\t{val}\n")

out.close()
print("wrote biosample_attributes_long.tsv")
PY



python3 - <<'PY'
import pandas as pd

df = pd.read_csv(
    "biosample_attributes_long.tsv",
    sep="\t",
    names=["BioSample","attr","val"],
    dtype=str
).dropna()

df["attr"] = df["attr"].str.strip()
df["val"]  = df["val"].str.strip()

keep = {
  "collection_date","geo_loc_name","geographic location",
  "lat_lon","latitude","longitude",
  "country","region","state","city",
  "host","isolation_source","strain","sample_name"
}

mask = df["attr"].isin(keep) | df["attr"].str.contains(
    r"(date|geo|locat|lat|lon|country|region|state|city)",
    case=False, na=False
)

wide = (
    df[mask]
    .groupby(["BioSample","attr"])["val"]
    .apply(lambda s: " | ".join(sorted(set(s))))
    .unstack(fill_value="")
)

wide.reset_index().to_csv("biosample_attributes_wide.tsv", sep="\t", index=False)
print("wrote biosample_attributes_wide.tsv")
PY


python3 - <<'PY'
import pandas as pd

runinfo = pd.read_csv("runinfo.tsv", sep="\t", dtype=str).fillna("")
bswide  = pd.read_csv("biosample_attributes_wide.tsv", sep="\t", dtype=str).fillna("")

m = runinfo.merge(bswide, on="BioSample", how="left")

front = [c for c in [
    "Run","Experiment","Sample","BioSample","BioProject",
    "ScientificName","TaxID",
    "collection_date","geo_loc_name","lat_lon"
] if c in m.columns]

m = m[front + [c for c in m.columns if c not in front]]
m.to_csv("sra_merged.tsv", sep="\t", index=False)

print("wrote sra_merged.tsv")
PY
