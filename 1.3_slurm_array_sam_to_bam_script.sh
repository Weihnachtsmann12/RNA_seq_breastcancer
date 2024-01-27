#!/bin/bash
#SBATCH --array=1-12
#SBATCH --time=13:00:00
#SBATCH --mem=40G
#SBATCH --cpus-per-task=5
#SBATCH --job-name=slurm_array_samtobam
#SBATCH --output=array_%J.out
#SBATCH --error=array_%J.err
#SBATCH --partition=pall

# define variables
WORKDIR="/data/users/mrubin/rnaseq_course"
OUTDIR="$WORKDIR/results"
SAMPLELIST="$WORKDIR/scripts/samplelist.tsv"

SAMPLE=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $1; exit}' $SAMPLELIST`
SAM=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $4; exit}' $SAMPLELIST`

OUTFILE="$OUTDIR/${SAMPLE}_sorted.bam"

############################

module load UHTS/Analysis/samtools/1.10
srun samtools view -hbS $SAM | samtools sort -o $OUTFILE
