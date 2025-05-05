

#conda create -n earlgrey_env python=3.9
conda activate earlgrey_env


#conda install -c bioconda -c conda-forge repeatmasker=4.1.8
#CONDA_VERBOSE=3 conda install -c bioconda -c conda-forge earlgrey



/home/rjp5nc/miniconda3/envs/earlgrey_env/share/earlgrey-4.2.4-0/earlGrey -g genome.fasta -s Daphniaobtusa -o earlgrey_output -l consensi.fa.classified


#consensi.fa.classified from repeatmodeller
#genome.fasta is the assembled fasta

