#!/bin/bash
set -euo pipefail

##########################################################################################################
# 07_bracken_abundance_cmds.sh
# Part of bioinformatics pipeline used in Eriksson et al (2025)
# Generate Bracken abundance estimation commands per sample
##########################################################################################################

kraken_report_dir="/path/to/kraken2_output"
kraken_db="/path/to/food_chl_lib" #kraken database
read_length=150
tax_level="S"
min_reads=10
output_cmd_file="bracken_cmds.txt"

> "$output_cmd_file"

for report in "$kraken_report_dir"/report.chl.*; do
    sample=$(basename "$report" | cut -d. -f3)
    output="${sample}.bracken"

    echo bracken -d $kraken_db -i $report -o $output -r $read_length -l $tax_level >> "$output_cmd_file"
done

echo "Bracken commands written to $output_cmd_file"
