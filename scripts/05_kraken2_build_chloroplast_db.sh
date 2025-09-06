#!/bin/bash

set -euo pipefail

##########################################################################################################
# 05_kraken2_build_chloroplast_db.sh
# Part of bioinformatics pipeline used in Eriksson et al (2025)
# Build chloroplast genome Kraken2 database for dietary classification
##########################################################################################################

#SBATCH --job-name=kraken2_build_chl
#SBATCH --partition=your_queue_name  # replace with your HPC queue/partition name
#SBATCH --time=0-23:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4

# Define database name
DB_NAME="food_chl_lib"

# Step 1: Download taxonomy
kraken2-build --download-taxonomy --db "$DB_NAME"

# Step 2: Add chloroplast genome FASTA files
for file in chloroplast_genomes/*.fasta; do
    kraken2-build --add-to-library "$file" --db "$DB_NAME"
done

# Step 3: Build the database
kraken2-build --build --db "$DB_NAME" --threads 4

echo "Chloroplast Kraken2 database build complete"
