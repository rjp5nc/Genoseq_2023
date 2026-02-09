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
conda init bash
source ~/.bashrc  # or restart your shell
conda activate pixy
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH

#conda install --yes -c conda-forge pixy
#conda install -c conda-forge pixy=2.0.0.beta11

#conda install -c bioconda htslib
#conda install -c bioconda vcftools
#conda install -c bioconda bcftools

cd /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv

file /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/trimmed10bp_allsites_usobtusa_renamed_annotated.vcf.gz


awk -F',' 'NR>1 {print $1 "\t" $6}' /project/berglandlab/Robert/UKSequencing2022_2024/old_stuff/2022_2024seqmetadata20250131.csv > pops.txt

awk 'NF==2' pops.txt > pops_filtered.txt

SEED=123  # change or remove for non-reproducible randomness

gawk -v SEED="$SEED" -v POPKEEP="P66" '
BEGIN{
  FS="[,\t ]+"
  OFS="\t"
  srand(SEED+0)
}
# -----------------------------
# 1) Read mitotype map: sample -> mito (A/B/C)
# -----------------------------
FNR==NR {
  if (FNR==1) next
  mito[$1]=$2
  next
}

# -----------------------------
# 2) Read genomic_types_v2.csv: sample -> superclone
#    Expected columns: ..., CloneA, Group
#    i.e. sample in $3, superclone in $4
# -----------------------------
ARGIND==2 {
  if (FNR==1) next
  s=$3
  g=$4
  gsub(/^"|"$/, "", s)
  gsub(/^"|"$/, "", g)
  if (s!="") super[s]=g
  next
}

# -----------------------------
# 3) Read pops_filtered.txt: keep only POPKEEP, gather subset
# -----------------------------
ARGIND==3 {
  s=$1
  pop=$2
  if (pop!=POPKEEP) next
  if (!(s in mito)) next
  if (!(s in super)) next

  # record sample list (for output)
  n++
  samp[n]=s

  # clone/superclone for this sample
  c=super[s]
  clone_of_sample[s]=c

  # original mitotype (for target counts)
  m=mito[s]
  orig_mito[s]=m
  target[m]++

  # clone sizes (in this filtered subset)
  clone_size[c]++
  next
}

END{
  # Ensure labels exist
  labs[1]="A"; labs[2]="B"; labs[3]="C"
  for (i=1;i<=3;i++) if (!(labs[i] in target)) target[labs[i]]=0

  # remaining counts we try to match across *samples*
  rem["A"]=target["A"]; rem["B"]=target["B"]; rem["C"]=target["C"]

  # Build clone list
  k=0
  for (c in clone_size) { k++; clones[k]=c }

  # Sort clones by decreasing size (simple O(k^2), fine for typical sizes)
  for (i=1;i<=k;i++){
    for (j=i+1;j<=k;j++){
      if (clone_size[clones[j]] > clone_size[clones[i]]) {
        tmp=clones[i]; clones[i]=clones[j]; clones[j]=tmp
      }
    }
  }

  # Assign each clone a mitotype label:
  # pick label that can "afford" size (rem>=size) with biggest remainder;
  # if none can afford, pick label with largest remainder (closest match).
  for (i=1;i<=k;i++){
    c=clones[i]
    s=clone_size[c]

    best=""
    bestScore=-1e18
    foundNonNeg=0

    # First pass: try only rem>=s
    for (li=1; li<=3; li++){
      L=labs[li]
      if (rem[L] >= s) {
        foundNonNeg=1
      }
    }

    for (li=1; li<=3; li++){
      L=labs[li]
      # score encourages using the label with the most "room" left
      # add tiny random jitter to break ties
      jitter = rand()*1e-6
      if (foundNonNeg) {
        if (rem[L] >= s) score = rem[L] + jitter
        else continue
      } else {
        score = rem[L] + jitter
      }
      if (score > bestScore) { bestScore=score; best=L }
    }

    clone_label[c]=best
    rem[best] -= s
  }

  # Output in same format as before: sample  pop_mito
  for (i=1;i<=n;i++){
    s=samp[i]
    c=clone_of_sample[s]
    L=clone_label[c]
    print s, POPKEEP "_" L
  }

  # stderr summary
  actualA=actualB=actualC=0
  for (c in clone_label){
    L=clone_label[c]
    if (L=="A") actualA += clone_size[c]
    if (L=="B") actualB += clone_size[c]
    if (L=="C") actualC += clone_size[c]
  }
  print "Target counts (samples): A="target["A"]" B="target["B"]" C="target["C"] > "/dev/stderr"
  print "Actual  counts (samples): A="actualA" B="actualB" C="actualC > "/dev/stderr"
}
' \
/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/allsites_mito/mito_types.csv \
/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/genomic_types_v2.csv \
pops_filtered.txt \
> pops_filtered_P66_withmitotype_rand_superclone.txt


join -t $'\t' \
  <(awk '$2 ~ /^P66_/' pops_filtered_withmitotype.txt | sort -k1,1) \
  <(sort -k1,1 pops_filtered_P66_withmitotype_rand_superclone.txt) \
| column -t



$CONDA_PREFIX/bin/bcftools query -l  /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/trimmed10bp_allsites_usobtusa_renamed_annotated.vcf.gz > vcf_samples.txt

grep -Ff vcf_samples.txt pops_filtered_P66_withmitotype_rand_superclone.txt > pops_both.txt

awk '$2 !~ /^PBO/' pops_filtered_P66_withmitotype_rand_superclone.txt > pops_fixed_withmitotype.noPBO.txt

grep -Ff <($CONDA_PREFIX/bin/bcftools query -l /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/trimmed10bp_allsites_usobtusa_renamed_annotated.vcf.gz) pops_fixed_withmitotype.noPBO.txt > pops_fixed_withmitotype2.txt


VCF="/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/trimmed10bp_allsites_usobtusa_renamed_annotated.vcf.gz"
POPS="/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/pops_fixed_withmitotype2.txt"

awk 'BEGIN{OFS="\t"} {print $1,$2}' pops_fixed_withmitotype2.txt \
  > pops_fixed_withmitotype.pixy.txt


# --- Output directory ---
OUTDIR="/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/results_pixy10000_withmitotype_random_superclone"
mkdir -p ${OUTDIR}

# --- Parameters ---
WINDOW=10000   # window size in bp, adjust as needed

chmod u+rwx /scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/results_pixy10000_withmitotype

