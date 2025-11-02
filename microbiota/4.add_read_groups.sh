#!/usr/bin/env bash
#SBATCH -J DownloadMap    # Job name
#SBATCH --ntasks=1        # Single task per job
#SBATCH --cpus-per-task=10 # Number of CPU cores per task
#SBATCH -N 1              # Run on one node
#SBATCH -t 1-20:00        # 10 hours runtime
#SBATCH --mem=100G        # Memory per node
#SBATCH -o /scratch/rjp5nc/outputerrors/down.%A_%a.out  # Standard output
#SBATCH -e /scratch/rjp5nc/outputerrors/down.%A_%a.err  # Standard error
#SBATCH -p standard       # Partition
#SBATCH --account=berglandlab

# Load modules
module load picard

# Move to directory
cd /scratch/rjp5nc/microbiota/chlorella/mapped_bam

# Force add read groups
java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
-I chlorella_ours_finalmap.bam \
-O /scratch/rjp5nc/microbiota/chlorella/mapped_bam/cholerlla_ours_finalmap_RG.bam \
-LB "library" \
-PL "ILLumina" \
-PU "platunit" \
-SM chlorella_ours 

# Index Bam files
java -jar $EBROOTPICARD/picard.jar BuildBamIndex \
-I /scratch/rjp5nc/microbiota/chlorella/mapped_bam/cholerlla_ours_finalmap_RG.bam