
src_dir="/project/berglandlab/daphnia_genus/short_read/murray_data/mito/new_mito.fasta"
dst_dir="/scratch/rjp5nc/UK2022_2024/consensusmito"

mkdir -p "$dst_dir"  # Create destination if it doesn't exist

for file in "$src_dir"/*; do
  base_file=$(basename "$file")
  dst_file="$dst_dir/$base_file"

  if [[ ! -e "$dst_file" ]]; then
    cp "$file" "$dst_file"
    echo "Copied: $base_file"
  else
    echo "Skipped (already exists): $base_file"
  fi
done
