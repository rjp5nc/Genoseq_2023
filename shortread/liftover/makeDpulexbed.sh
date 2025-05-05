awk '
    /^>/ {
        if (seqname != "") {
            for (i = 0; i < length(seq); i++) {
                base = substr(seq, i+1, 1);
                print seqname "\t" i "\t" i+1 "\t" base;
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
        for (i = 0; i < length(seq); i++) {
            base = substr(seq, i+1, 1);
            print seqname "\t" i "\t" i+1 "\t" base;
        }
    }
' /scratch/rjp5nc/Reference_genomes/post_kraken/eu_pulex_totalHiCwithallbestgapclosed.clean.fa > /scratch/rjp5nc/Reference_genomes/eu_pulex_totalHiCwithallbestgapclosed.clean.allbases.bed


head -n 20 /scratch/rjp5nc/Reference_genomes/eu_pulex_totalHiCwithallbestgapclosed.clean.allbases.bed

tail -n 20 /scratch/rjp5nc/Reference_genomes/eu_pulex_totalHiCwithallbestgapclosed.clean.allbases.bed







awk '
    /^>/ {
        if (seqname != "") {
            for (i = 0; i < length(seq); i++) {
                base = substr(seq, i+1, 1);
                print seqname "\t" i "\t" i+1 "\t" base;
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
        for (i = 0; i < length(seq); i++) {
            base = substr(seq, i+1, 1);
            print seqname "\t" i "\t" i+1 "\t" base;
        }
    }
' /scratch/rjp5nc/Reference_genomes/post_kraken/us_pulex_ref_kap4.fa > /scratch/rjp5nc/Reference_genomes/us_pulex_ref_kap4.allbases.bed

head -n 20 /scratch/rjp5nc/Reference_genomes/us_pulex_ref_kap4.allbases.bed

awk '{print $1, $(NF-2), $(NF-1), $NF}' OFS='\t' /scratch/rjp5nc/Reference_genomes/us_pulex_ref_kap4.allbases.bed > /scratch/rjp5nc/Reference_genomes/us_pulex_ref_kap4.allbases.cleaned.bed

head -n 20 /scratch/rjp5nc/Reference_genomes/us_pulex_ref_kap4.allbases.cleaned.bed


awk '
    /^>/ {
        if (seqname != "") {
            for (i = 0; i < length(seq); i++) {
                base = substr(seq, i+1, 1);
                print seqname "\t" i "\t" i+1 "\t" base;
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
        for (i = 0; i < length(seq); i++) {
            base = substr(seq, i+1, 1);
            print seqname "\t" i "\t" i+1 "\t" base;
        }
    }
' /scratch/rjp5nc/Reference_genomes/post_kraken/assembly.hap2_onlydaps.fasta > /scratch/rjp5nc/Reference_genomes/assembly.hap2_onlydaps.allbases.bed

head -n 20 /scratch/rjp5nc/Reference_genomes/assembly.hap2_onlydaps.allbases.bed



awk '
    /^>/ {
        if (seqname != "") {
            for (i = 0; i < length(seq); i++) {
                base = substr(seq, i+1, 1);
                print seqname "\t" i "\t" i+1 "\t" base;
            }
        }
        seqname = substr($1, 2);
        seq = "";
        next;
    }
    {
        seq = seq $0;
    }
    END {
        for (i = 0; i < length(seq); i++) {
            base = substr(seq, i+1, 1);
            print seqname "\t" i "\t" i+1 "\t" base;
        }
    }
' /scratch/rjp5nc/Reference_genomes/post_kraken/Daphnia_ambigua_Q001_genome.fa > /scratch/rjp5nc/Reference_genomes/Daphnia_ambigua_Q001_genome.allbases.bed

head -n 20 /scratch/rjp5nc/Reference_genomes/Daphnia_ambigua_Q001_genome.allbases.bed


awk '
    /^>/ {
        if (seqname != "") {
            for (i = 0; i < length(seq); i++) {
                base = substr(seq, i+1, 1);
                print seqname "\t" i "\t" i+1 "\t" base;
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
        for (i = 0; i < length(seq); i++) {
            base = substr(seq, i+1, 1);
            print seqname "\t" i "\t" i+1 "\t" base;
        }
    }
' /scratch/rjp5nc/Reference_genomes/post_kraken/US_obtusa_onlydaps.fa > /scratch/rjp5nc/Reference_genomes/US_obtusa_onlydaps.allbases.bed


head -n 20 /scratch/rjp5nc/Reference_genomes/US_obtusa_onlydaps.allbases.bed

awk '{print $1, $(NF-2), $(NF-1), $NF}' OFS='\t' /scratch/rjp5nc/Reference_genomes/US_obtusa_onlydaps.allbases.bed > /scratch/rjp5nc/Reference_genomes/US_obtusa_onlydaps.allbases.cleaned.bed

head -n 20 /scratch/rjp5nc/Reference_genomes/US_obtusa_onlydaps.allbases.cleaned.bed




cd /scratch/rjp5nc/Reference_genomes


head -n 20 /scratch/rjp5nc/Reference_genomes/Daphnia_ambigua_Q001_genome.allbases.bed



input=us_pulex_ref_kap4.allbases.cleaned
input=US_obtusa_onlydaps.allbases.cleaned
input=Daphnia_ambigua_Q001_genome.allbases

input=assembly.hap2_onlydaps.allbases

awk '{print $1, $2, $3, $4, $1","$2","$4}' OFS='\t' $input.bed > $input.slimmed.bed

head -n 20 $input.slimmed.bed






awk '{print $1, $2, $3, $4, $1","$2","$4}' OFS='\t' /scratch/rjp5nc/Reference_genomes/eu_pulex_totalHiCwithallbestgapclosed.clean.allbases.bed > /scratch/rjp5nc/Reference_genomes/eu_pulex_totalHiCwithallbestgapclosed.clean.allbases.slimmed.bed
