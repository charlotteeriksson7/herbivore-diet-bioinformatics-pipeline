#!/bin/bash

set -euo pipefail

##########################################################################################################
# 09_samtools_pipeline.sh
# Part of bioinformatics pipeline used in Eriksson et al (2025)
# Full pipeline — Convert SAM to sorted BAM, generate stats and MAPQ distribution,
# filter reads, and run flagstat pre- and post-filtering
##########################################################################################################

# SLURM settings
#SBATCH --job-name=samtools_pipeline
#SBATCH --time=1-00:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --partition=your_queue_name  # replace with your HPC queue/partition name
#SBATCH --nodes=1

# Settings
MEMORY="4G"
THREADS=4
MAPQ_THRESHOLD=48
INPUT_DIR="data/"
OUTPUT_BASE_DIR="samtools_results"
LOG_DIR="logs"

# Output subdirectories
OUTPUT_DIR="${OUTPUT_BASE_DIR}/mapping"
STATS_DIR="${OUTPUT_BASE_DIR}/stats"
MAPQ_DIR="${OUTPUT_BASE_DIR}/mapq"
FLAGSTAT_DIR="${OUTPUT_BASE_DIR}/flagstat"
FILTERED_DIR="${OUTPUT_BASE_DIR}/filtered"
FLAGSTAT_FILTERED_DIR="${OUTPUT_BASE_DIR}/flagstat_filtered"


mkdir -p "$OUTPUT_DIR" "$STATS_DIR" "$MAPQ_DIR" "$FLAGSTAT_DIR" "$FILTERED_DIR" "$FLAGSTAT_FILTERED_DIR" "$LOG_DIR"

echo "Starting samtools filtering pipeline..."
echo "SLURM Job ID: $SLURM_JOB_ID"
echo "Input directory: $INPUT_DIR"

for SAM in "$INPUT_DIR"/*_min.sam; do
    BASENAME=$(basename "$SAM" .sam)

# Main loop
    SORTED_BAM="${OUTPUT_DIR}/${BASENAME}_sorted.bam"
    STATS_FILE="${STATS_DIR}/${BASENAME}_sumstats.txt"
    MAPQ_FILE="${MAPQ_DIR}/${BASENAME}_mapq.txt"
    FLAGSTAT_FILE="${FLAGSTAT_DIR}/${BASENAME}_flagstat.txt"
    FILTERED_BAM="${FILTERED_DIR}/${BASENAME}_sorted_filt.bam"
    FLAGSTAT_FILTERED_FILE="${FLAGSTAT_FILTERED_DIR}/${BASENAME}_filt_flagstat.txt"

    echo "Processing $BASENAME"

    # Step 1: Convert SAM to sorted BAM
    samtools view -bS "$SAM" | \
        samtools sort -@ "$THREADS" -m "$MEMORY" -o "$SORTED_BAM" -

    # Step 2: Summary statistics
    samtools stats "$SORTED_BAM" | grep ^SN | cut -f 2- > "$STATS_FILE"

    # Step 3: MAPQ distribution
    samtools view "$SORTED_BAM" | cut -f 5 | sort | uniq -c | sort -n | \
        awk '{printf("MAPQ:%s\t%d\n",$2,$1);}' > "$MAPQ_FILE"

    # Step 4: Flagstat pre-filtering
    samtools flagstat "$SORTED_BAM" > "$FLAGSTAT_FILE"

    # Step 5: Filter reads based on flags and MAPQ
    samtools view -b -F 2308 -f 0x2 -q "$MAPQ_THRESHOLD" "$SORTED_BAM" > "$FILTERED_BAM"

    # Step 6: Flagstat post-filtering
    samtools flagstat "$FILTERED_BAM" > "$FLAGSTAT_FILTERED_FILE"

    echo "Finished $BASENAME"
    echo "  → Sorted BAM: $SORTED_BAM"
    echo "  → Stats: $STATS_FILE"
    echo "  → MAPQ: $MAPQ_FILE"
    echo "  → Pre-filter flagstat: $FLAGSTAT_FILE"
    echo "  → Filtered BAM: $FILTERED_BAM"
    echo "  → Post-filter flagstat: $FLAGSTAT_FILTERED_FILE"
done


echo "All samples processed successfully."
