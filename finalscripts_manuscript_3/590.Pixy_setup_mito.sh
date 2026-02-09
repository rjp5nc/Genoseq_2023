#!/usr/bin/env bash
#
#SBATCH -J pixy # A single job name for the array
#SBATCH --ntasks-per-node=20 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-92:00  ### 48 hours
#SBATCH --mem 80G
#SBATCH -o /scratch/rjp5nc/err/pixy.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/err/pixy.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab


#ijob -A berglandlab -c10 -p standard --mem=40G

#conda create -n pixy -c conda-forge -c bioconda pixy -y
# activate the environment
conda activate pixy

#conda install --yes -c conda-forge pixy
#conda install -c conda-forge pixy=2.0.0.beta11

#conda install -c bioconda htslib
#conda install -c bioconda vcftools


module load bcftools
cd /scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/pixy_all_wholemito




meta=/project/berglandlab/Robert/UKSequencing2022_2024/old_stuff/metadata_with_clone.csv
samps=/scratch/rjp5nc/UK2022_2024/daphnia_phylo/Sample_ID_Species_merged_20251227.csv
out=/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/pixy_all_wholemito/clone_accuratelocation.csv

awk -F, -v OFS=',' '
  function deq(x){ gsub(/^"|"$/, "", x); return x }

  BEGIN{ print "Well","clone","accuratelocation" }

  # ---- 1) metadata_with_clone.csv: keep all rows ----
  FNR==1 && NR==1 { next }   # skip meta header
  NR==FNR {
    well  = deq($2)
    clone = deq($5)
    acc   = deq($9)

    # drop blank clones if you want (optional)
    # if (tolower(clone) ~ /^blank$/) next

    print well, clone, acc
    seen[well]=1
    next
  }

  # ---- 2) species file: add rows for Sample_ID_old prefixes ----
  FNR==1 { next }            # skip species header
  {
    sid  = deq($3)           # Sample_ID (SRR...)
    sold = deq($4)           # Sample_ID_old (EBG-62 etc)

    if (sold ~ /^(FS|RAP|bdw|PYR|EBG)/) {
      pref = sold
      sub(/[^A-Za-z].*$/, "", pref)  # keep just leading letters (EBG, PYR, bdw, etc)

      # if this SRR already exists in metadata, still output an override row;
      # downstream you can prefer these rows by last-one-wins if needed.
      print sid, sold, pref
    }
  }
' "$meta" "$samps" > "$out"

echo "Wrote: $out"
tail "$out"
wc -l "$out"



# cp /scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/gatk_gvcf/usdobtusa_mito_joint_Full.vcf.gz /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/



awk -F',' 'NR>1 {print $1 "\t" $3}' /scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/pixy_all_wholemito/clone_accuratelocation.csv > pops.txt

awk 'NF==2' pops.txt > pops_filtered.txt



bcftools query -l  /scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/gatk_gvcf/usdobtusa_mito_joint_Full.vcf.gz > vcf_samples.txt





awk -F'[,\t ]+' '
NR==FNR {
  if (FNR>1) mito[$1]=$2; next
}
($1 in mito) {
  print $1, $2 "_" mito[$1]
}
' /scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/mito_types.csv pops_filtered.txt > pops_filtered_withmitotype.txt




$CONDA_PREFIX/bin/bcftools query -l  /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/usdobtusa_mito_joint_Full.vcf.gz > vcf_samples.txt

grep -Ff vcf_samples.txt pops_filtered_withmitotype.txt > pops_both.txt


grep -Ff <(bcftools query -l /scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/gatk_gvcf/usdobtusa_mito_joint_Full.vcf.gz) pops_both.txt > pops_fixed.txt


VCF="/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/gatk_gvcf/usdobtusa_mito_joint_Full.vcf.gz"
POPS="/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/pixy_all_wholemito/pops_fixed.txt"

# --- Output directory ---
OUTDIR="/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/pixy_all_wholemito/"
mkdir -p ${OUTDIR}

# --- Parameters ---

cd /scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/pixy_all_wholemito/
tabix -p vcf /scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/gatk_gvcf/usdobtusa_mito_joint_Full.vcf.gz
chmod u+rwx /scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/pixy_whole_mitoresults/

awk '{print $1}' pops_fixed.txt > pops_samples.txt


bcftools query -l /scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/gatk_gvcf/usdobtusa_mito_joint_Full.vcf.gz > vcf_samples.txt

# Keep only matching samples
grep -F -x -f vcf_samples.txt pops_samples.txt > pops_samples_filtered.txt

bcftools view -S pops_samples_filtered.txt -Oz -o usdobtusa_mito_genotyped_subset.vcf.gz /scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/gatk_gvcf/usdobtusa_mito_joint_Full.vcf.gz

grep -F -w -f pops_samples_filtered.txt pops_fixed.txt > pops_fixed_filtered.txt

tabix -p vcf usdobtusa_mito_genotyped_subset.vcf.gz

bcftools query -f '[%DP\t]\n' usdobtusa_mito_genotyped_subset.vcf.gz | \
awk '
{for(i=1;i<=NF;i++){sum[i]+=$i; count[i]++}}
END{
    for(i=1;i<=NF;i++){
        mean=sum[i]/count[i];
        print i, mean
    }
}' > sample_mean_depth.txt

bcftools query -l usdobtusa_mito_genotyped_subset.vcf.gz > sample_names.txt
paste sample_names.txt sample_mean_depth.txt > sample_depths_mapped.txt

awk '$3<30 {print $1}' sample_depths_mapped.txt > low_coverage_samples.txt

bcftools view -s ^$(paste -sd, low_coverage_samples.txt) \
    -Oz -o usdobtusa_mito_filtered.vcf.gz --force-samples \
    usdobtusa_mito_genotyped_subset.vcf.gz

bcftools index -t usdobtusa_mito_filtered.vcf.gz
bcftools query -l usdobtusa_mito_filtered.vcf.gz > filtered_samples.txt
grep -Ff filtered_samples.txt pops_fixed_filtered.txt > pops_filtered_for_pixy.txt

bcftools query -l usdobtusa_mito_filtered.vcf.gz > vcf_samples.txt
awk 'NR==FNR {a[$1]=1; next} ($1 in a) {print $1 "\t" $2}' \
  vcf_samples.txt pops_fixed_filtered.txt > pops_pixy_ready.txt


bcftools view -H usdobtusa_mito_filtered.vcf.gz | head -n 20 | awk '{print $1,$2,$4,$5}' | column -t


bcftools view \
  -r CM028013.1:1421-2129 \
  -Oz -o cox1.vcf.gz \
  usdobtusa_mito_filtered.vcf.gz

bcftools index -t cox1.vcf.gz


bcftools query -f'[%GT\n]' cox1.vcf.gz | head


pixy --stats pi fst dxy \
--vcf cox1.vcf.gz \
--populations pops_pixy_ready.txt \
--window_size 2129 \
--n_cores 40 \
--output_folder /scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/pixy_whole_mitoresults/ \
--bypass_invariant_check \
--output_prefix pixy_COI


pixy --stats pi fst dxy \
--vcf usdobtusa_mito_filtered.vcf.gz \
--populations pops_pixy_ready.txt \
--window_size 14685 \
--n_cores 40 \
--output_folder /scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/pixy_whole_mitoresults/ \
--bypass_invariant_check \
--output_prefix pixy

# --- Done ---
echo "✅ Pixy run finished. Outputs in ${OUTDIR}/"


#awk '{print $1 > $2".txt"}' pops_fixed_filtered.txt

#for popfile in P58.txt P62.txt P63.txt P66.txt Gilmer.txt; do
#    popname=$(basename "$popfile" .txt)
#    vcftools --gzvcf usdobtusa_mito_genotyped_subset.vcf.gz \
#             --keep "$popfile" \
#             --site-pi \
#             --out "${popname}_pi"
#done