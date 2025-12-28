#!/usr/bin/env bash
#SBATCH -J DownloadMap    # Job name
#SBATCH --ntasks=1        # Single task per job
#SBATCH --cpus-per-task=10 # Number of CPU cores per task
#SBATCH -N 1              # Run on one node
#SBATCH -t 1-20:00        # 10 hours runtime
#SBATCH --mem=25G        # Memory per node
#SBATCH -o /scratch/rjp5nc/outputerrors/down.%A_%a.out  # Standard output
#SBATCH -e /scratch/rjp5nc/outputerrors/down.%A_%a.err  # Standard error
#SBATCH -p standard       # Partition
#SBATCH --account=berglandlab
#SBATCH --array=22-286

# Load modules
module load picard

#find /scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/bam/ -name "*_finalmap.bam" | sort > /scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/bam_files_newref.txt

# Parameters
parameterFile="/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/bam_files_newref.txt"
wd="/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/final_mitobam_rg3"

# Extract sample name
samp=`sed -n ${SLURM_ARRAY_TASK_ID}p $parameterFile`
out2=`echo $samp | sed 's/_finalmap.bam//'`
out=$(echo "$out2" | sed 's#/scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/bam/##')

echo "Adding read groups -" "Sample:" $SLURM_ARRAY_TASK_ID

# Move to directory
cd /scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/final_mitobam_rg3/

# Force add read groups
java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
-I ${samp} \
-O $wd/${out}finalmap_RG.bam \
-LB "library" \
-PL "ILLumina" \
-PU "platunit" \
-SM ${out} 

# Index Bam files
java -jar $EBROOTPICARD/picard.jar BuildBamIndex \
-I $wd/${out}finalmap_RG.bam



#cp /scratch/rjp5nc/UK2022_2024/final_mitobam_rg3/*.b /scratch/rjp5nc/UK2022_2024/NA1_Dobtusa/final_mitobam_rg3/

