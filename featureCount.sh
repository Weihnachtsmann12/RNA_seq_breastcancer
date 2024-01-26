#!/bin/bash

#SBATCH --job-name=featureCounts_paired_S1
#SBATCH --output=featureCounts_paired_S1.out
#SBATCH --error=featureCounts_paired_S1.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --time=15:00:00
#SBATCH --mem=10G
#SBATCH --partition=pall


bam_path="/data/users/mrubin/rnaseq_course/results/"
ref_path="/data/users/mrubin/rnaseq_course/reference_genome/"

module load UHTS/Analysis/subread/2.0.1
featureCounts -t exon -g gene_id -p -S 2 -a ${ref_path}Homo_sapiens.GRCh38.110.gtf -o featurecounts_paired_S1.txt ${bam_path}*sorted.bam