# Repeat Sequence Analysis Pipeline
A comprehensive Python pipeline for repeat sequence analysis using BLAST, with quality filtering, statistical analysis, and visualization capabilities.

## Overview
This pipeline identifies and characterizes repetitive sequences in genomic data through BLAST-based homology search. It features:
- Random sequence extraction from reference genomes (when no query sequence is provided)
- Quality filtering of repeat sequences with customizable thresholds
- Statistical analysis of repeat units (copy number, identity, length)
- High-quality repeat sequence identification
- Visualization of repeat distribution across chromosomes
- Automatic result organization and logging
- Resumable analysis (skip completed steps)

## Table of Contents
1. [Installation](#installation)
2. [Dependencies](#dependencies)
3. [Configuration](#configuration)
4. [Usage](#usage)
5. [Output Files](#output-files)
6. [Parameters](#parameters)
7. [Troubleshooting](#troubleshooting)

## Installation
```bash
# Clone or download the script
git clone <repository-url>
cd repeat-analysis

# Make the script executable
chmod +x repeat_analysis.py

# Install Python dependencies
pip install -r requirements.txt
```

## Dependencies
### System Dependencies
| Tool | Purpose | Minimum Version |
|------|---------|-----------------|
| blastn | Nucleotide BLAST search | 2.10+ |
| makeblastdb | BLAST database creation | 2.10+ |
| blastdbcmd | BLAST database verification | 2.10+ |
| samtools | Sequence extraction (optional) | 1.9+ |
| seqkit | Sequence extraction (optional) | 2.0+ |
| awk | Text processing | GNU Awk 5.0+ |
| sort | File sorting | GNU coreutils 8.0+ |

### Python Dependencies
See `requirements.txt`:
```txt
pandas>=1.3.0
numpy>=1.21.0
matplotlib>=3.4.0
```

## Configuration
### Create Configuration Template
```bash
python repeat_analysis.py --create-config
```
This generates `repeat_analysis_config.ini` with default parameters.

### Configuration Parameters
Edit the INI file to set your analysis parameters:
```ini
[analysis]
# Input query sequence file (optional, random extraction if empty)
query = 

# Reference genome FASTA (required if query is empty)
reference_genome = 

# BLAST database name/path (required)
genome_db = 

# Output file prefix (optional)
output_prefix = 

# Computational resources
threads = 4

# BLAST filtering
e_value = 1e-50
min_identity = 80
min_length = 100

# Repeat unit filtering
min_copies = 1000
min_avg_identity = 90
min_avg_length = 200

# Output settings
top_n = 10

# Visualization
chr_lengths = 

# Random sequence extraction (when query is empty)
extract_length = 5000000
extract_regions = 1
```

## Usage
### Basic Usage
```bash
# Run with default config file
python repeat_analysis.py

# Run with custom config file
python repeat_analysis.py my_config.ini

# Force re-analysis (ignore existing results)
python repeat_analysis.py my_config.ini --force

# Create config template only
python repeat_analysis.py --create-config
```

### Typical Workflow
1. Create configuration file with `--create-config`
2. Edit configuration to set your input files and parameters
3. Run the script with your config file
4. Check output directory for results and visualizations

## Output Files
The pipeline creates a structured output directory with the following subdirectories:

| Directory | Content |
|-----------|---------|
| `blast_results/` | Raw and sorted BLAST output files |
| `repeat_analysis/` | Repeat unit statistics and filtered results |
| `visualization/` | Chromosome distribution plots (if enabled) |
| `fasta_files/` | Extracted random query sequences (if used) |
| `config_files/` | Generated chromosome length files and config backups |

### Key Output Files
- `<prefix>_blast.txt`: Raw BLAST results
- `<prefix>_repeat_units.txt`: All identified repeat units
- `<prefix>_high_quality_repeats.txt`: Filtered high-quality repeats
- `<prefix>_top_repeats.txt`: Top N high-quality repeats by copy number
- `<prefix>_statistics_report.txt`: Summary statistics
- `repeat_analysis.log`: Detailed analysis log

## Parameters
### Command Line Arguments
| Argument | Description |
|----------|-------------|
| `config_file` | Path to configuration file (default: repeat_analysis_config.ini) |
| `--create-config` | Generate configuration template and exit |
| `--force` | Force re-run of all analysis steps |

### Filtering Parameters
| Parameter | Description | Default |
|-----------|-------------|---------|
| `e_value` | BLAST e-value threshold | 1e-50 |
| `min_identity` | Minimum sequence identity (%) | 80 |
| `min_length` | Minimum alignment length (bp) | 100 |
| `min_copies` | Minimum copy number for repeat units | 1000 |
| `min_avg_identity` | Minimum average identity for high-quality repeats (%) | 90 |
| `min_avg_length` | Minimum average length for high-quality repeats (bp) | 200 |
| `top_n` | Number of top repeats to report | 10 |

## Troubleshooting
### Common Issues
1. **BLAST database errors**:
   - Ensure `genome_db` points to valid BLAST database files or a FASTA file
   - The script will attempt to create a BLAST database if a FASTA file is provided
   
2. **Missing dependencies**:
   - Install missing system tools with your package manager (e.g., `apt install ncbi-blast+ samtools`)
   - Install Python dependencies with `pip install -r requirements.txt`
   
3. **No repeat units found**:
   - Check filtering parameters (reduce `min_copies` or `min_identity`)
   - Verify BLAST database and query sequences are compatible
   - Check log file for detailed error messages

4. **Memory issues with large datasets**:
   - Increase available memory or reduce `extract_length`
   - Process data in smaller chunks or use fewer threads

### Log File
Check `repeat_analysis.log` for detailed error messages and processing information.

## Notes
- The pipeline supports both standard INI config files and legacy key-value format files
- Results are resumable - the script skips completed steps unless `--force` is used
- Chromosome length files are automatically generated if not provided
- Random sequence extraction requires a reference genome FASTA file
- High-quality repeats are defined by both average identity and length thresholds

## License
[Add your license information here]

## Contact
[Add contact information or issue tracker link here]