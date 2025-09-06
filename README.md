## Herbivore Diet Bioinformatics Pipeline

This repository contains all bioinformatics scripts used in the metagenomics and targeted capture portion of this study:

Comparing accuracy and biases of DNA metabarcoding, hybridization capture, and metagenomic sequencing for quantifying herbivore diets  

> *Currently under peer review*  
> Author: Charlotte Eriksson · Contact: charlotte.eriksson@oregonstate.edu


## Repository Structure

```herbivore-diet-bioinformatics-pipeline/
├── README.md # This file
├── LICENSE # GNU General Public License v3.0
  └── scripts/ # Analysis scripts
    ├── 01_trim_and_clean.sh # fastp: Adapter/quality trimming
    ├── 02_qc_fastqc.sh # FastQC: Quality assessment 
    ├── 03_kraken2_build_reference_db.sh # Kraken2: Full reference database
    ├── 04_kraken2_classify.sh # Kraken2: Taxonomic classification
    ├── 05_kraken2_build_chloroplast_db.sh # Kraken2: Chloroplast database
    ├── 06_kraken2_chl_classify.sh # Kraken2: Chloroplast classification
    ├── 07_bracken_abundance.sh # Bracken: Relative abundance estimation
    ├── 08_minimap2_mapping.sh # Minimap2: Chloroplast read mapping
    ├── 09_samtools_pipeline.sh # Samtools: Filtering mapped reads
    ├── 09b_samtools_dedup_pipeline.sh # Samtools: Filtering with deduplication
    └── 10_genome_coverage_analysis.sh # BEDTools: Genome coverage calculations

## Pipeline Summary

### 1. Trimming and Cleaning
- **Tool**: `fastp`
- **Script**: `01_trim_and_clean.sh`
- **Purpose**: Trim adapters, remove low-quality reads (Phred < 30), and filter reads < 50 bp.

### 2. Quality Control
- **Tool**: `FastQC`
- **Script**: `02_fastqc.sh`
- **Purpose**: Generate per-sample quality reports post-trimming.

### 3. Taxonomic Classification (Full Database)
- **Tool**: `Kraken2`
- **Scripts**:
  - `03_kraken2_build_reference_db.sh`: Builds custom database with RefSeq, nt, UniVec, and deer genomes.
  - `04_kraken2_classify.sh`: Classifies reads using the full database (`--paired` mode).

### 4. Chloroplast Classification
- **Tool**: `Kraken2`
- **Scripts**:
  - `05_kraken2_build_chloroplast_db.sh`: Builds database with chloroplast genomes of 24 plant species.
  - `06_kraken2_chl_classify.sh`: Classifies using chloroplast database with confidence thresholds:
    - `0.5` for metagenomic samples
    - `0.05` for hybrid capture samples

### 5. Abundance Estimation
- **Tool**: `Bracken`
- **Script**: `07_bracken_abundance.sh`: Estimates relative abundance
- **Filters**:
  - Discards taxa <0.1% of total reads (metagenomic)
  - Discards taxa <1% of total reads (hybrid capture)

### 6. Read Mapping
- **Tool**: `Minimap2`
- **Script**: `08_minimap2_mapping.sh`
- **Purpose**: Aligns reads to full chloroplast genomes.

### 7. Read Filtering
- **Tool**: `Samtools`
- **Scripts**:
  - `09_samtools_pipeline.sh`: Filters out unmapped, secondary, and supplementary reads (`-F 2308`)
  - `09b_samtools_dedup_pipeline.sh`: Also deduplicates reads (for hybrid capture)
- **Additional Filters**: MAPQ ≥ 48 and properly paired reads (`-q 48 -f 0x2`)

### 8. Genome Coverage Calculation
- **Tool**: `BEDTools`
- **Script**: `10_genome_coverage_analysis.sh`
- **Purpose**: Calculates per-base chloroplast genome coverage, normalized by genome length

---

## Software Dependencies

Install using conda:

```bash
conda create -n diet-pipeline fastp fastqc kraken2 bracken minimap2 samtools bedtools -c bioconda

## Input Data Notes

Raw FASTQ files and reference genomes are not included due to size constraints.

Chloroplast reference genome accessions are listed in Supplemental Table S1 of the manuscript.

Processed data will be uploaded to Dryad upon publication.

## License

This project is licensed under the GNU General Public License v3.0.
See the LICENSE file for full details

## Contact

For questions, feedback, or contributions, please contact:

Charlotte Eriksson
charlotte.eriksson@oregonstate.edu
