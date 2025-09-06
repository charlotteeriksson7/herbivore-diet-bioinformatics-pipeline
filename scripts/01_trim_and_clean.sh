#!/bin/bash

set -euo pipefail

##########################################################################################################
# 01_trim_and_clean.sh
# Part of bioinformatics pipeline used in Eriksson et al (2025)
# Trim adapters, remove low quality reads (Phred < 30) and filter out sequences < 50 bp using fastp
##########################################################################################################

#SBATCH --partition=your_queue_name # replace with your HPC queue/partition name
#SBATCH --job-name=fastp_trim
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

# Create output directory for cleaned files
mkdir -p cleaned

# Loop over all R1 fastq files and run fastp paired with corresponding R2
for R1 in *R1*.fastq; do
  R2=${R1/R1/R2}  # Assumes R1 and R2 have matching names
  OUT1=cleaned/${R1%.fastq}.cleaned.fastq
  OUT2=cleaned/${R2%.fastq}.cleaned.fastq

  echo "Processing $R1 and $R2 ..."
  fastp -i "$R1" -I "$R2" -o "$OUT1" -O "$OUT2" -q 30 -l 50 -w 3 --detect_adapter_for_pe --trim_poly_g
done
