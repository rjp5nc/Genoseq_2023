#!/usr/bin/env bash
#
#SBATCH -J merge # A single job name for the array
#SBATCH --ntasks-per-node=2 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 2-10:00 # 10 hours
#SBATCH --mem 20G
#SBATCH -o /scratch/rjp5nc/erroroutputs/down.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/erroroutputs/down.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab


# Configuration
csv_file_old="/scratch/rjp5nc/UK2022_2024/daphnia_phylo/uspulexsampsmove.txt"
csv_file="/scratch/rjp5nc/UK2022_2024/daphnia_phylo/uspulexsampsmove_amended.txt"
target_dir="/scratch/rjp5nc/UK2022_2024/daphnia_phylo/gvcf/uspulex_chr"
column_name="Library_Name"
# Get the column index for the sample_id

sed 's/\r/\n/g' "$csv_file_old" > "$csv_file"

cd /project/berglandlab/connor/mapping_pulex_nam/mapping/gvcf/pulex

col_num=$(awk -F',' -v col="$column_name" '
NR==1 {
  for (i=1; i<=NF; i++) {
    gsub(/^"|"$/, "", $i)
    if ($i == col) {
      print i
      exit
    }
  }
}' "$csv_file")

if [ -z "$col_num" ]; then
  echo "Column '$column_name' not found."
  exit 1
fi

# Loop over each sample ID
tail -n +2 "$csv_file" | awk -F',' -v col="$col_num" '{gsub(/^"|"$/, "", $col); print $col}' | while read -r sample_id; do
  find . -type f -name "*$sample_id*" | while read -r file; do
    rel_path="${file#./}"  # remove leading ./ from path
    dest_path="$target_dir/$rel_path"
    mkdir -p "$(dirname "$dest_path")"
    cp "$file" "$dest_path"
    echo "Copied: $file -> $dest_path"
  done
done