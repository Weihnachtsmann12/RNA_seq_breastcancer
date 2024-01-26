#!/usr/bin/env bash
#SBATCH --job-name=hisat2_index
#SBATCH --output=hisat2_index.out
#SBATCH --error=hisat2_index.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=3:00:00
#SBATCH --partition=pall

module load UHTS/Aligner/hisat/2.2.1
srun hisat2-build -p 16 Homo_sapiens.GRCh38.dna.primary_assembly.fa Homo_sapiens.GRCh38.dna.primary_assembly.fa