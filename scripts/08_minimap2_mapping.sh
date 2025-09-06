#!/bin/bash
set -euo pipefail

##########################################################################################################
# 08_minimap2_mapping.sh
# Part of bioinformatics pipeline used in Eriksson et al (2025)
# Generate Minimap2 mapping commands
##########################################################################################################

input_dir="/path/to/fastq_files"
reference_fasta="/path/to/chloroplast_library.fasta"
output_cmd_file="minimap2_cmds.txt"

for r1 in "$input_dir"/*_R1.fastq; do
    sample=$(basename "$r1" _R1.fastq)
    r2="${input_dir}/${sample}_R2.fastq"

    echo "minimap2 -ax sr $reference_fasta $r1 $r2 > ${sample}_min.sam" >> "$output_cmd_file"
done

echo "Commands written to $output_cmd_file"
