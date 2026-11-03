module load miniforge
conda create -y -n mito -c bioconda -c conda-forge python=3.8 diamond
conda activate mito
pip install tiara
tiara --help | head -3
diamond version

mkdir -p /scratch/rjp5nc/mitoDB && cd /scratch/rjp5nc/mitoDB
wget https://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/mitochondrion.1.protein.faa.gz
diamond makedb --in mitochondrion.1.protein.faa.gz -d refseq_mito