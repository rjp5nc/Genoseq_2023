

conda create -n repeatmodeler_env -c bioconda -c conda-forge repeatmodeler perl=5.22.0
conda activate repeatmodeler_new

seqtk seq dap.contigs.fasta > cleaned_dap.contigs.fasta

makeblastdb -in cleaned_dap.contigs.fasta -out my_db -dbtype nucl -title my_db -parse_seqids
BuildDatabase -name my_db dap.contigs.fasta
RepeatModeler -database my_db -pa 4