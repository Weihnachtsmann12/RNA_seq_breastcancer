#!/bin/bash
#SBATCH --array=1-12
#SBATCH --time=13:00:00
#SBATCH --mem=10G
#SBATCH --cpus-per-task=4
#SBATCH --job-name=slurm_array_mapping
#SBATCH --output=array_%J.out
#SBATCH --error=array_%J.err
#SBATCH --partition=pall

# define variables
WORKDIR="/data/users/mrubin/rnaseq_course"
OUTDIR="$WORKDIR/results"
SAMPLELIST="$WORKDIR/scripts/samplelist.tsv"

SAMPLE=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $1; exit}' $SAMPLELIST`
READ1=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $2; exit}' $SAMPLELIST`
READ2=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $3; exit}' $SAMPLELIST`

OUTFILE="$OUTDIR/${SAMPLE}.sam"

############################


mkdir -p $OUTDIR

module load UHTS/Aligner/hisat/2.2.1
srun hisat2 -x /data/users/mrubin/rnaseq_course/reference_genome/hisat2/Homo_sapiens.GRCh38.dna.primary_assembly.fa -1 $READ1 -2 $READ2 -S $OUTFILE