#!/usr/bin/env bash
#
#SBATCH -J run_modeler # A single job name for the array
#SBATCH --ntasks-per-node=15 # multi core
#SBATCH -N 1 # on one node
#SBATCH -t 5-00:00 # 1 days
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

wd="/scratch/rjp5nc/removedups/eu_dobtusa/"
cd ${wd}

#seqtk seq assembly.hap2_onlydaps.fasta > cleaned_dap.contigs.fasta
#makeblastdb -in cleaned_dap.contigs.fasta -out my_db -dbtype nucl -title my_db -parse_seqids
#BuildDatabase -name my_db assembly.hap2_onlydaps.fasta

#RepeatModeler -database my_db -pa 15


RepeatModeler -database my_db -pa 15 -recoverDir RM_776028.SunMar231202032025
