#!/usr/bin/env bash
#
#SBATCH -J run_blast # A single job name for the array
#SBATCH --ntasks-per-node=10 # multi core
#SBATCH -N 1 # on one node
#SBATCH -t 0-03:00 # 3 hours
#SBATCH --mem 40G
#SBATCH -o /scratch/rjp5nc/Canu_error/blast%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/Canu_error/blast%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

#conda create -n mitohifi_env python=3.8
conda activate mitohifi_env

module load hifiasm samtools blast


#git clone https://github.com/RemiAllio/MitoFinder.git
#cd /home/rjp5nc/bin/MitoFinder
#python3 setup.py install

#export PATH="/home/rjp5nc/miniconda3/envs/mitohifi_env/bin:$PATH"
#export PATH="/home/rjp5nc/bin/MitoHiFi/src:$PATH"

#seqtk seq -F '#' m84128_250121_222443_s2.hifi_reads.bc2104.fq.gz > m84128_250121_222443_s2.hifi_reads.bc2104.fastq

python3 /home/rjp5nc/bin/MitoHiFi/src/mitohifi.py -r none -f /scratch/rjp5nc/HMW/HMWDNAElvis3/m84128_250121_222443_s2.hifi_reads.bc2104.fastq -g /scratch/rjp5nc/HMW/HMWDNAElvis3/us_obtusa.gb -t 20 -o /scratch/rjp5nc/HMW/HMWDNAElvis3/mitohifi_output

python3 /home/rjp5nc/bin/MitoHiFi/src/mitohifi.py -r /scratch/rjp5nc/HMW/HMWDNAElvis3/m84128_250121_222443_s2.hifi_reads.bc2104.fastq -f /scratch/rjp5nc/HMW/HMWDNAElvis3/us_obtusa_mito.fa -t 20 -g /scratch/rjp5nc/HMW/HMWDNAElvis3/us_mito.gb -o /scratch/rjp5nc/HMW/HMWDNAElvis3/mitohifi_output




makeblastdb -in assembly.p_ctg_forblast.fa -out obtusadbhap1 -dbtype nucl -title obtusadbhap1 -parse_seqids
blastn -db obtusadbhap1 -query us_mito.fasta -outfmt "6 qacc sacc length pident qstart qend sstart send qseq" -evalue 1e-5 -out obtusamito.txt



#h1tg000076l hap1 mitochondrial 

awk '/^>h1tg000076l$/ {flag=1; print; next} /^>/ {flag=0} flag {print}' assembly.p_ctg_forblast.fa > h1tg000076lhap1.fasta


awk '/^>/ {if (seq) print seq; print; seq=""; next} {seq = seq $0} END {print seq}' h1tg000076lhap1.fasta | \
  awk 'BEGIN {start=10720; end=25388} {if (NR >= start && NR <= end) print substr($0, 1)}' > output_file.fasta