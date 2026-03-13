
conda create -y -n pysra -c bioconda -c conda-forge pysradb
conda activate pysra

cd /scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/


pysradb srr-to-srp $(cat sra_USobtusa.txt) \
  | awk 'NR>1 {print $2}' \
  | sort -u > srp_list.txt

> sra_metadata_all.tsv
first=1

while read -r SRP; do
  [ -z "$SRP" ] && continue
  if [ $first -eq 1 ]; then
    pysradb metadata --desc --expand --detailed "$SRP" > sra_metadata_all.tsv
    first=0
  else
    pysradb metadata --desc --expand --detailed "$SRP" | tail -n +2 >> sra_metadata_all.tsv
  fi
done < srp_list.txt

head -n 1 sra_metadata_all.tsv | tr '\t' '\n' | grep -i -E 'geo|loc|country|state|province|site'


awk 'NR==FNR {keep[$1]=1; next}
     FNR==1 || ($1 in keep)' \
    sra_USobtusa.txt sra_metadata_all.tsv \
    > sra_metadata_filtered.tsv

cut -f1 sra_metadata_filtered.tsv | tail -n +2 | sort > present.txt
sort sra_USobtusa.txt > expected.txt
comm -23 expected.txt present.txt

sed 's/\t/,/g' sra_metadata_filtered.tsv > sra_metadata_filtered.csv

