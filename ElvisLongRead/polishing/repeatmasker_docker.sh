# ijob -A berglandlab -c4 -p standard --mem=10G

### RepeatMakser image
### https://github.com/Dfam-consortium/TETools

module load apptainer

### install docker
singularity pull dfam-tetools-latest.sif docker://dfam/tetools:latest
curl -sSLO https://github.com/Dfam-consortium/TETools/raw/master/dfam-tetools.sh
chmod +x dfam-tetools.sh
/scratch/rjp5nc/dfamtools/dfam-tetools.sh -- RepeatModeler -h
/scratch/rjp5nc/dfamtools/dfam-tetools.sh -- RepeatMasker -h

### make library


#finished us obtusa, eu obtusa, us ambigua, eu pulex


ref=/scratch/rjp5nc/Reference_genomes/post_kraken/US_obtusa_onlydaps.fa
classified=/scratch/rjp5nc/removedups/us_dobtusa/us_dobtusa/RM_653473.SunMar301757312025/consensi.fa.classified

cd /scratch/rjp5nc/dfamtools/

### try running RepeatMasker
~/dfam-tetools.sh -- RepeatMasker -h
/scratch/rjp5nc/dfamtools/dfam-tetools.sh \
-- RepeatMasker \
-e rmblast \
-pa 4 \
-qq \
-lib $classified \
-dir /scratch/rjp5nc/removedups/us_dobtusa/ \
-html \
$ref