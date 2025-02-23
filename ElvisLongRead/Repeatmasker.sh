#!/usr/bin/env bash
#
#SBATCH -J run_modeler # A single job name for the array
#SBATCH --ntasks-per-node=20 # multi core
#SBATCH -N 1 # on one node
#SBATCH -t 1-00:00 # 1 days
#SBATCH --mem 100G
#SBATCH -o /scratch/rjp5nc/outputerrors/m%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/outputerrors/m%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab





#conda create -n quast_env python=3.9
#conda activate quast_env
#conda install -c bioconda quast

#quast.py /scratch/rjp5nc/assemblecontigs/dap.contigs.fasta -o /scratch/rjp5nc/assemblecontigs


#conda create -n repeatmodeler_env -c bioconda -c conda-forge repeatmodeler perl=5.22.0
conda activate repeatmodeler_new
#seqtk seq dap.contigs.fasta > cleaned_dap.contigs.fasta
#makeblastdb -in cleaned_dap.contigs.fasta -out my_db -dbtype nucl -title my_db -parse_seqids
#BuildDatabase -name my_db dap.contigs.fasta

wd="/scratch/rjp5nc/removedups"
cd ${wd}

RepeatModeler -database my_db -pa 20
