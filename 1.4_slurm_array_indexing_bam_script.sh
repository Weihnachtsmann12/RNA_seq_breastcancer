#!/bin/bash
#SBATCH --array=1-12
#SBATCH --time=10:00:00
#SBATCH --mem=5G
#SBATCH --cpus-per-task=1
#SBATCH --job-name=slurm_array_indexbam
#SBATCH --output=array_%J.out
#SBATCH --error=array_%J.err
#SBATCH --partition=pall

# define variables
WORKDIR="/data/users/mrubin/rnaseq_course"
SAMPLELIST="$WORKDIR/scripts/samplelist.tsv"

BAM=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $5; exit}' $SAMPLELIST`

############################

module load UHTS/Analysis/samtools/1.10
srun samtools index $BAM
