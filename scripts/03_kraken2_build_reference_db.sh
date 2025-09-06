#!/bin/bash

set -euo pipefail

##########################################################################################################
# 03_kraken2_build_reference_db.sh
# Part of bioinformatics pipeline used in Eriksson et al (2025)
# Build reference Kraken2 database for assessing proportion of plant, deer, and microbe DNA
##########################################################################################################

#SBATCH --job-name=kraken2_build
#SBATCH --partition=your_queue_name  # replace with your HPC queue/partition name
#SBATCH --time=2-00:00:00
#SBATCH --mem=80G
#SBATCH --cpus-per-task=5

# Define database name
DB_NAME="reference_db"

# Create directory
mkdir -p $DB_NAME

# Download RefSeq libraries
kraken2-build --download-library archaea --db $DB_NAME
kraken2-build --download-library bacteria --db $DB_NAME
kraken2-build --download-library fungi --db $DB_NAME
kraken2-build --download-library human --db $DB_NAME
kraken2-build --download-library plant --db $DB_NAME
kraken2-build --download-library plasmid --db $DB_NAME
kraken2-build --download-library protozoa --db $DB_NAME
kraken2-build --download-library UniVec_Core --db $DB_NAME
kraken2-build --download-library viral --db $DB_NAME

# Download nt (NCBI nucleotide database)
kraken2-build --download-library nt --db $DB_NAME

# Add custom genomes
kraken2-build --add-to-library genomes/GCA_020976825.1_BYU_Ohem_2021.fna --db $DB_NAME   # Mule deer
kraken2-build --add-to-library genomes/GCA_002102435.1_Ovir.te_1.0.fna --db $DB_NAME     # White-tailed deer

# Build the full database
kraken2-build --build --db $DB_NAME --threads 5
