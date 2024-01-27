#!/usr/bin/env bash
#SBATCH --job-name=fastqc_samples
#SBATCH --output=fastqc_samples.out
#SBATCH --error=fastqc_samples.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=0:45:00
#SBATCH --partition=pall

module add UHTS/Quality_control/fastqc/0.11.9
srun fastqc --extract /data/users/mrubin/rnaseq_course/fastqc/fastq_files/* --threads 1