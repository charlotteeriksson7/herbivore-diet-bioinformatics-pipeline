#!/bin/bash

set -euo pipefail

################################################################################
# 10_genome_coverage_analysis.sh
# Part of the bioinformatics pipeline (Eriksson et al. 2025)
# Calculates genome coverage statistics per sample from filtered BAM files
################################################################################

# SLURM directives (if submitting with sbatch)
#SBATCH --job-name=genome_coverage
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --partition=your_queue_name  # replace with your HPC queue/partition name
#SBATCH --nodes=1

# Settings
INPUT_DIR="samtools_results/filtered"
OUTPUT_DIR="samtools_results/genome_coverage"
THREADS=2

# Create output directory
mkdir -p "$OUTPUT_DIR"

echo "Starting genome coverage analysis..."

# Loop over all filtered BAM files
for BAM in "$INPUT_DIR"/*_sorted_filt.bam; do
    BASENAME=$(basename "$BAM" _sorted_filt.bam)
    echo "Processing $BASENAME"

    # Step 1: Index BAM
    samtools index "$BAM"

    # Step 2: Calculate genome coverage at each base
    COVERAGE_FILE="${OUTPUT_DIR}/${BASENAME}_genomecov.txt"
    bedtools genomecov -d -ibam "$BAM" > "$COVERAGE_FILE"

    # Step 3: Count bases with >0 coverage per contig
    PERBASE_COUNT_FILE="${OUTPUT_DIR}/${BASENAME}.perbase.count.txt"
    awk -F"\t" '$3>0{print $1}' "$COVERAGE_FILE" | sort | uniq -c > "$PERBASE_COUNT_FILE"

    # Step 4: Extract contig lengths from BAM header
    LENGTH_FILE="${OUTPUT_DIR}/${BASENAME}.lengths.genome"
    samtools view -H "$BAM" | \
        awk '/^@SQ/ {for (i=1; i<=NF; i++) if ($i ~ /^SN:/) sn=substr($i,4); else if ($i ~ /^LN:/) ln=substr($i,4);} {print sn "\t" ln}' > "$LENGTH_FILE"

    # Step 5: Calculate genome coverage proportion per contig
    PROPORTION_FILE="${OUTPUT_DIR}/${BASENAME}.genomeproportion"
    while IFS=$'\t' read -r contig length; do
        covered=$(awk -v c="$contig" '$2==c {print $1}' "$PERBASE_COUNT_FILE")
        covered=${covered:-0}  # Default to 0 if not found
        proportion=$(echo "scale=5; $covered / $length" | bc)
        echo "${contig},${proportion}"
    done < "$LENGTH_FILE" > "$PROPORTION_FILE"

    echo "  Coverage file:         $COVERAGE_FILE"
    echo "  Covered bases count:   $PERBASE_COUNT_FILE"
    echo "  Genome lengths:        $LENGTH_FILE"
    echo "  Genome proportions:    $PROPORTION_FILE"
done

echo "Genome coverage analysis complete for all samples."
