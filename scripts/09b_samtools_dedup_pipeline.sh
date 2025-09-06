#!/bin/bash

set -euo pipefail

##########################################################################################################
# 09b_samtools_dedup_pipeline.sh
# Part of bioinformatics pipeline used in Eriksson et al (2025)
# Deduplication version — for hybridization capture data only
# Converts SAM to BAM, deduplicates reads using fixmate + markdup,
# filters by flags and MAPQ, and generates summary stats
##########################################################################################################

# SLURM settings
#SBATCH --job-name=samtools_dedup
#SBATCH --time=1-00:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --partition=your_queue_name  # replace with your HPC queue/partition name
#SBATCH --nodes=1

# Settings
MEMORY="4G"
THREADS=4
MAPQ_THRESHOLD=48
INPUT_DIR="data/"
OUTPUT_BASE_DIR="samtools_results"

# Output subdirectories
OUTPUT_DIR="${OUTPUT_BASE_DIR}/mapping_dedup"
STATS_DIR="${OUTPUT_BASE_DIR}/stats_dedup"
MAPQ_DIR="${OUTPUT_BASE_DIR}/mapq_dedup"
FLAGSTAT_DIR="${OUTPUT_BASE_DIR}/flagstat_dedup"
FILTERED_DIR="${OUTPUT_BASE_DIR}/filtered_dedup"
FLAGSTAT_FILTERED_DIR="${OUTPUT_BASE_DIR}/flagstat_filtered_dedup"

# Create output directories
mkdir -p "$OUTPUT_DIR" "$STATS_DIR" "$MAPQ_DIR" "$FLAGSTAT_DIR" "$FILTERED_DIR" "$FLAGSTAT_FILTERED_DIR"

echo "Starting SAMtools deduplication + filtering pipeline..."
echo "SLURM Job ID: $SLURM_JOB_ID"
echo "Input directory: $INPUT_DIR"

for SAM in "$INPUT_DIR"/*_min.sam; do
    BASENAME=$(basename "$SAM" .sam)

    echo "Processing $BASENAME"

    # Define filenames
    SORTED_BAM="${OUTPUT_DIR}/${BASENAME}_sorted.bam"
    NAME_SORTED_BAM="${OUTPUT_DIR}/${BASENAME}_name_sorted.bam"
    FIXMATE_BAM="${OUTPUT_DIR}/${BASENAME}_fixmate.bam"
    COORD_SORTED_BAM="${OUTPUT_DIR}/${BASENAME}_coord_sorted.bam"
    DEDUP_BAM="${OUTPUT_DIR}/${BASENAME}_dedup.bam"
    FILTERED_BAM="${FILTERED_DIR}/${BASENAME}_sorted_filt.bam"

    STATS_FILE="${STATS_DIR}/${BASENAME}_sumstats.txt"
    MAPQ_FILE="${MAPQ_DIR}/${BASENAME}_mapq.txt"
    FLAGSTAT_FILE="${FLAGSTAT_DIR}/${BASENAME}_flagstat.txt"
    FLAGSTAT_FILTERED_FILE="${FLAGSTAT_FILTERED_DIR}/${BASENAME}_filt_flagstat.txt"

    # Step 1: Convert SAM to coordinate-sorted BAM
    samtools view -@ "$THREADS" -bS "$SAM" | \
        samtools sort -@ "$THREADS" -m "$MEMORY" -o "$SORTED_BAM" -

    # Step 2: Name sort
    samtools sort -n -@ "$THREADS" -o "$NAME_SORTED_BAM" "$SORTED_BAM"

    # Step 3: Fixmate
    samtools fixmate -m "$NAME_SORTED_BAM" "$FIXMATE_BAM"

    # Step 4: Coordinate sort again
    samtools sort -@ "$THREADS" -o "$COORD_SORTED_BAM" "$FIXMATE_BAM"

    # Step 5: Remove PCR duplicates
    samtools markdup -r "$COORD_SORTED_BAM" "$DEDUP_BAM"

    # Step 6: Summary statistics on deduplicated BAM
    samtools stats "$DEDUP_BAM" | grep ^SN | cut -f 2- > "$STATS_FILE"

    # Step 7: MAPQ distribution
    samtools view "$DEDUP_BAM" | cut -f 5 | sort | uniq -c | sort -n | \
        awk '{printf("MAPQ:%s\t%d\n",$2,$1);}' > "$MAPQ_FILE"

    # Step 8: Flagstat pre-filtering
    samtools flagstat "$DEDUP_BAM" > "$FLAGSTAT_FILE"

    # Step 9: Apply filtering (MAPQ ≥ 48, properly paired, remove unmapped, secondary, supplementary)
    samtools view -b -F 2308 -f 0x2 -q "$MAPQ_THRESHOLD" "$DEDUP_BAM" > "$FILTERED_BAM"

    # Step 10: Flagstat post-filtering
    samtools flagstat "$FILTERED_BAM" > "$FLAGSTAT_FILTERED_FILE"

    echo "Finished $BASENAME"
    echo "  Deduplicated BAM: $DEDUP_BAM"
    echo "  Filtered BAM: $FILTERED_BAM"
done

echo "All samples processed with deduplication successfully."
