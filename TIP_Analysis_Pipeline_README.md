# TIP Analysis Pipeline
A comprehensive bash pipeline for analyzing TIP (Transposable Element Insertion Polymorphism) data, including quality control, genome alignment, read filtering, tag analysis, and SNP calling for male/female sample comparisons.

## Overview
This pipeline processes paired-end sequencing data to identify sex-specific TIPs and SNPs, with a focus on transposable element analysis. It combines bioinformatics tools for data processing and custom Python/awk scripts for statistical analysis and filtering.

## Pipeline Workflow
The pipeline executes the following steps sequentially:
1. **Data Quality Control** (fastp) - Adapter trimming and quality filtering
2. **Genome Alignment** (bwa mem) - Mapping reads to reference genome
3. **BAM Filtering** (bamtools/samtools) - Retaining uniquely mapped reads
4. **Read Position Extraction** - Capturing genomic coordinates of filtered reads
5. **Tag Cluster Analysis** - Identifying sex-specific/shared tag clusters
6. **Specific Tag Filtering** - Finding tags containing transposon terminal sequences
7. **BAM Integrity Check** - Validating BAM files before variant calling
8. **SNP Calling** (bcftools) - Variant detection from aligned reads
9. **SNP Classification** - Identifying sex-specific and shared polymorphisms
10. **SNP Analysis** - Calculating MAF and statistical summaries

## Dependencies
### Required Software
| Tool | Version Requirement | Purpose |
|------|---------------------|---------|
| bash | ≥4.0 | Pipeline execution |
| fastp | ≥0.23.0 | Read quality control |
| bwa | ≥0.7.17 | Genome alignment |
| samtools | ≥1.15 | BAM file processing/indexing |
| bamtools | ≥2.5.1 | BAM filtering |
| bcftools | ≥1.15 | SNP calling and filtering |
| python3 | ≥3.6 | Tag cluster analysis |
| awk | GNU awk (gawk) ≥5.0 | Text processing |
| coreutils | Latest | Basic file operations (wc, sort, etc.) |

### Python Libraries
No external Python libraries required (uses only standard libraries: `csv`, `collections`).

### Configuration Requirements
A `config_1.ini` file with the following parameters must be present in the working directory:
```ini
# Sample information
MALE_SAMPLE=male_sample_name
FEMALE_SAMPLE=female_sample_name

# Directory paths
WORKDIR=/path/to/output_directory
RAW_DATA_DIR=/path/to/raw_fastq_files
REF_GENOME=/path/to/reference_genome.fasta

# Quality thresholds
MIN_LENGTH=30
THREADS=16
MAPQ_THRESHOLD=30
MAX_MISMATCH=2
CLUSTER_DISTANCE=50
SPECIFIC_SEQUENCE=TGCAACCTGTGTTCCGGTGCC

# SNP calling parameters
SNP_BASEQ=20
SNP_MAPQ=30
MAX_DEPTH=1000
PLOIDY=2
SNP_QUALITY=20
MIN_DP=5
MAF_THRESHOLD=0.05

# Output file names
MALE_SPECIFIC_FILE=male_specific.tsv
FEMALE_SPECIFIC_FILE=female_specific.tsv
```

## Installation
1. Install required tools via conda (recommended):
   ```bash
   conda create -n tip_analysis -c bioconda fastp bwa samtools bamtools bcftools python=3.8 gawk
   conda activate tip_analysis
   ```

2. Clone/download this pipeline:
   ```bash
   git clone <repository-url>
   cd tip_analysis
   ```

3. Prepare your `config_1.ini` file with appropriate paths and parameters

## Usage
### Basic Execution
```bash
# Make the script executable
chmod +x TIP_analysis.sh

# Run the pipeline
./TIP_analysis.sh
```

### Expected Input
- Paired-end FASTQ files (gzipped) for male and female samples:
  - `${MALE_SAMPLE}_R1.fq.gz` and `${MALE_SAMPLE}_R2.fq.gz`
  - `${FEMALE_SAMPLE}_R1.fq.gz` and `${FEMALE_SAMPLE}_R2.fq.gz`
- Reference genome FASTA file (indexed with `bwa index`)
- Configuration file (`config_1.ini`) with all parameters defined

### Output Structure
The pipeline creates the following directory structure in `$WORKDIR`:
```
$WORKDIR/
├── clean_data/          # Quality-filtered FASTQ + fastp reports
├── alignment/
│   ├── raw/             # Unfiltered BAM files
│   └── unique/          # Filtered (unique mapping) BAM files
├── read_position/       # Read coordinate files + tag cluster results
├── specific_tags/       # Transposon sequence-containing tags
├── vcf/
│   ├── raw/             # Unfiltered VCF files
│   └── filtered/        # Quality-filtered VCF files
├── snp/                 # SNP classification and MAF analysis results
└── logs/                # Log files for each processing step
```

## Output File Descriptions
| File Path | Description |
|-----------|-------------|
| `read_position/tag_clusters.tsv` | All tag clusters with sex classification |
| `read_position/*_specific_tags.tsv` | Sex-specific tag clusters |
| `specific_tags/*_with_sequence.tsv` | Reads containing transposon terminal sequences |
| `specific_tags/*_tag_stats.tsv` | Statistics of matched reads per tag |
| `snp/both_polymorphic.tsv` | Shared SNPs with MAF values for both sexes |
| `snp/*_only_polymorphic.tsv` | Sex-specific SNPs with MAF calculations |
| `snp/*.stats.txt` | Detailed SNP statistics from bcftools |

## Troubleshooting
1. **BAM index errors**: Ensure samtools version ≥1.15 (uses CSI indexes by default)
2. **SNP calling failures**: Check BAM file integrity with `samtools quickcheck`
3. **Tag analysis issues**: Verify read position files exist and are properly formatted
4. **Memory issues**: Reduce `THREADS` parameter or increase system resources
5. **Missing sequence matches**: Adjust `SPECIFIC_SEQUENCE` or `CLUSTER_DISTANCE` parameters

## Notes
- All intermediate files are retained for debugging/downstream analysis
- Log files in `logs/` directory contain detailed error messages for troubleshooting
- Adjust parameters in `config_1.ini` based on your specific dataset characteristics
- For large datasets, ensure sufficient disk space (minimum 50GB recommended)
- The pipeline is designed for diploid organisms (adjust `PLOIDY` parameter if needed)

## License
This pipeline is provided as-is for research purposes. Please cite appropriate tools (fastp, BWA, samtools, bcftools) when using this pipeline in publications.

## Contact
For pipeline-related issues, please open an issue in the repository or contact the developer.