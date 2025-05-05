./liftOver input.bed hg19ToHg38.over.chain output.bed unMapped.bed

#samtools faidx eu_pulex_totalHiCwithallbestgapclosed.clean.fa
#cut -f1,2 eu_pulex_totalHiCwithallbestgapclosed.clean.fa.fai | awk '{print $1 "\t0\t" $2}' > eu_pulex_all_chr.bed


conda create -n liftover-env -c bioconda ucsc-liftover
conda activate liftover-env


bedtools makewindows -b /scratch/rjp5nc/Reference_genomes/assembly.hap2_onlydaps/assembly.hap2_onlydaps_all_scaffolds.bed -w 100000 > /scratch/rjp5nc/Reference_genomes/assembly.hap2_onlydaps/assembly.hap2_onlydaps_all_scaffolds_split100000.bed

# I need to redo liftover the other way around
#The map.chain file has the old genome as the target and the new genome
#as the query.

~/miniconda3/envs/liftover-env/bin/liftOver /scratch/rjp5nc/Reference_genomes/assembly.hap2_onlydaps/assembly.hap2_onlydaps_all_scaffolds_split100000.bed \
 /scratch/rjp5nc/lastz/eu_obtusa/chainnet/liftover.chain \
 /scratch/rjp5nc/outputtest/output100000.bed \
 /scratch/rjp5nc/outputtest/unmapped100000.bed


 less  /scratch/rjp5nc/lastz/eu_obtusa/chainnet/liftover.chain