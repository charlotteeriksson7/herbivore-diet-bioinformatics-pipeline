#!/bin/bash

set -euo pipefail

##########################################################################################################
# 02_qc_fastqc.sh
# Part of bioinformatics pipeline used in Eriksson et al (2025)
# Assessing quality with FastQC 
##########################################################################################################

#SBATCH --partition=your_queue_name # replace with your HPC queue/partition name
#SBATCH --job-name=fastqc
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=8G

# Create output folder
mkdir -p fastqc_results

# Run FastQC on all cleaned fastq files
fastqc cleaned/*.cleaned.fastq -o fastqc_results -t 6

# Run MultiQC to aggregate results
multiqc fastqc_results -o fastqc_results
