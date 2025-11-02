#!/usr/bin/env bash
#
#SBATCH -J windowfasta # A single job name for the array
#SBATCH --ntasks-per-node=5 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-2:00:00 # 8 hours
#SBATCH --mem 60G
#SBATCH -o /scratch/rjp5nc/err/gatk.chrom.%A_%a.out # Standard output
#SBATCH -e /scratch/rjp5nc/err/gatk.chrom.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab


# Input FASTA file

awk -v window=1000000 'BEGIN {
    print "id,seqnames,start,end,width";
}
$0 ~ /^>/ {
    if (seqname != "") {
        len = length(seq);
        for (start = 1; start <= len; start += window) {
            end = start + window - 1;
            if (end > len) end = len;
            width = end - start + 1;
            print id++ "," seqname "," start "," end "," width;
        }
    }
    seqname = substr($0, 2);
    seq = "";
    next;
}
{
    seq = seq $0;
}
END {
    if (seqname != "") {
        len = length(seq);
        for (start = 1; start <= len; start += window) {
            end = start + window - 1;
            if (end > len) end = len;
            width = end - start + 1;
            print id++ "," seqname "," start "," end "," width;
        }
    }
}' /scratch/rjp5nc/Reference_genomes/post_kraken/Daphnia_ambigua_Q001_genome.fa > /scratch/rjp5nc/UK2022_2024/daphnia_phylo/interval_DBI_paramList_ambigua.txt















awk -v window=1000000 '
BEGIN {
    print "id,seqnames,start,end,width";
    id = 1;
}
/^>/ {
    if (seqname != "" && seqlen > 0) {
        for (start = 1; start <= seqlen; start += window) {
            end = start + window - 1;
            if (end > seqlen) end = seqlen;
            width = end - start + 1;
            print id++ "," seqname "," start "," end "," width;
        }
    }
    # Extract only the contig ID from the header (up to first space or parenthesis)
    match($0, /^>([^ ]+)/, arr);
    seqname = arr[1];
    seqlen = 0;
    next;
}
{
    seqlen += length($0);
}
END {
    if (seqname != "" && seqlen > 0) {
        for (start = 1; start <= seqlen; start += window) {
            end = start + window - 1;
            if (end > seqlen) end = seqlen;
            width = end - start + 1;
            print id++ "," seqname "," start "," end "," width;
        }
    }
}
' /scratch/rjp5nc/Reference_genomes/post_kraken/Daphnia_ambigua_Q001_genome.fa > /scratch/rjp5nc/UK2022_2024/daphnia_phylo/interval_DBI_paramList_ambigua3.txt