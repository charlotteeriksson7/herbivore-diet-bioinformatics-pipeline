#!/bin/bash
set -euo pipefail
 
##########################################################################################################
# 04_kraken2_classify_cmds.sh
# Part of bioinformatics pipeline used in Eriksson et al (2025)
# Generate Kraken2 classification commands for all FASTQ files
##########################################################################################################


input_dir="/path/to/fastq"
kraken_db="path/to/reference_db"
threads=3
output_cmd_file="kraken2_cmds.txt"

for r1 in "$input_dir"/*_R1.fastq; do
    filename=$(basename "$r1")
    sample_prefix="${filename%_R1.fastq}"
    r2="$input_dir/${sample_prefix}_R2.fastq"

    echo "kraken2 --threads $threads --report report.chl.$sample_prefix --paired $r1 $r2 --db $kraken_db --output Run.kraken.out.chl.$sample_prefix --memory-mapping" >> "$output_cmd_file"
done
