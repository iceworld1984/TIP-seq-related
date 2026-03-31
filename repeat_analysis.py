#!/usr/bin/env python3
"""
Repeat Sequence Analysis Script - Python Version
Features:
1. Repeat sequence analysis (based on BLAST)
2. Support for random sequence extraction from reference genome as query sequences
3. Detailed logging output
4. High-quality repeat sequence filtering and visualization
5. Chromosome distribution map generation
6. Automatic result file organization
7. Support for resuming analysis from existing results
Note:
Normal first run: python script.py config.ini
Resume analysis (using existing results): python script.py config.ini
Force re-analysis: python script.py config.ini --force
"""

import os
import sys
import logging
import configparser
import argparse
import subprocess
import tempfile
import random
import shutil
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

class RepeatAnalyzer:
    """Repeat sequence analyzer"""
    
    def __init__(self):
        # Default parameters
        self.default_config = {
            'query': '',
            'reference_genome': '',
            'genome_db': '',
            'output_prefix': '',
            'threads': '4',
            'e_value': '1e-50',
            'min_identity': '80',
            'min_length': '100',
            'min_copies': '1000',
            'min_avg_identity': '90',
            'min_avg_length': '200',
            'top_n': '10',
            'chr_lengths': '',
            'extract_length': '5000000',  # Default extraction length: 5Mb
            'extract_regions': '1'        # Default number of regions: 1
        }
        
        self.config = self.default_config.copy()
        self.output_dir = None
        self.chrom_lengths = {}
        self.existing_results = {}  # Store existing result file paths
        self.setup_logging()
    
    def setup_logging(self):
        """Setup logging"""
        logging.basicConfig(
            level=logging.INFO,
            format='[%(levelname)s] %(asctime)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S',
            handlers=[
                logging.StreamHandler(),
                logging.FileHandler('repeat_analysis.log', mode='w', encoding='utf-8')
            ]
        )
        self.logger = logging.getLogger(__name__)
    
    def log_info(self, message: str):
        """Info log"""
        self.logger.info(message)
    
    def log_warn(self, message: str):
        """Warning log"""
        self.logger.warning(message)
    
    def log_error(self, message: str):
        """Error log"""
        self.logger.error(message)
        sys.exit(1)
    
    def create_output_directory(self):
        """Create output directory"""
        if self.config['query']:
            base_name = Path(self.config['query']).stem
        else:
            base_name = "random_queries"
        
        self.output_dir = Path(base_name)
        self.output_dir.mkdir(exist_ok=True)
        
        # Create result subdirectory structure
        subdirs = ['blast_results', 'repeat_analysis', 'visualization', 'fasta_files', 'config_files']
        for subdir in subdirs:
            (self.output_dir / subdir).mkdir(exist_ok=True)
        
        # Update output prefix to path within directory
        self.config['output_prefix'] = str(self.output_dir / base_name)
        self.log_info(f"Created output directory: {self.output_dir}")
        self.log_info(f"Output file prefix: {self.config['output_prefix']}")
    
    def get_output_path(self, filename: str, subdir: str = None) -> str:
        """Get output file path, optionally save to subdirectory"""
        if subdir:
            # Ensure subdirectory exists
            subdir_path = self.output_dir / subdir
            subdir_path.mkdir(exist_ok=True)
            return str(subdir_path / filename)
        else:
            # Return base path
            return f"{self.config['output_prefix']}_{filename}"


    def parse_arguments(self):
        """Parse command line arguments"""
        parser = argparse.ArgumentParser(description='Repeat sequence analysis script')
        parser.add_argument('config_file', nargs='?', default='repeat_analysis_config.ini',
                          help='Configuration file path (default: repeat_analysis_config.ini)')
        parser.add_argument('--create-config', action='store_true',
                          help='Create configuration file template and exit')
        parser.add_argument('--force', action='store_true',
                          help='Force re-run all steps, ignoring existing results')
        return parser.parse_args()
    
    def create_config_template(self):
        """Create configuration file template"""
        config_content = """# Repeat Sequence Analysis Configuration File
# Please modify the following parameters according to your needs

[analysis]
# Input query sequence file (optional, if empty, random sequences will be extracted from genome)
query = 

# Reference genome FASTA file (used for random sequence extraction, required if query is empty)
reference_genome = 

# BLAST database name (required, can be database name or FASTA file)
genome_db = 

# Output file prefix (optional, defaults to query sequence file name)
output_prefix = 

# Number of threads
threads = 4

# BLAST e-value threshold
e_value = 1e-50

# Minimum sequence identity (%)
min_identity = 80

# Minimum sequence length (bp)
min_length = 100

# Minimum copy number requirement
min_copies = 1000

# High-quality repeat sequence standard - average identity (%)
min_avg_identity = 90

# High-quality repeat sequence standard - average length (bp)
min_avg_length = 200

# Number of top repeat sequences to return
top_n = 10

# Chromosome length file (for plotting, optional, if not provided, will be auto-generated from reference genome)
chr_lengths = 

# Random extraction sequence length (bp, used when query is empty)
extract_length = 5000000

# Number of random extraction regions (used when query is empty)
extract_regions = 1
"""
        
        with open('repeat_analysis_config.ini', 'w', encoding='utf-8') as f:
            f.write(config_content)
        
        self.log_info("Configuration file template created: repeat_analysis_config.ini")
    
    def parse_legacy_config(self, config_file: str):
        """Parse legacy format configuration file (without section headers)"""
        self.log_info(f"Detected legacy format configuration file, attempting to parse: {config_file}")
        
        config_dict = {}
        try:
            with open(config_file, 'r', encoding='utf-8') as f:
                for line in f:
                    line = line.strip()
                    # Skip empty lines and comments
                    if not line or line.startswith('#') or line.startswith(';'):
                        continue
                    
                    # Parse key-value pairs
                    if '=' in line:
                        key, value = line.split('=', 1)
                        key = key.strip()
                        value = value.strip()
                        
                        # Remove possible quotes
                        if (value.startswith('"') and value.endswith('"')) or \
                           (value.startswith("'") and value.endswith("'")):
                            value = value[1:-1]
                        
                        config_dict[key] = value
            
            # Update configuration
            for key in self.default_config:
                if key in config_dict:
                    self.config[key] = config_dict[key]
            
            self.log_info("Successfully parsed legacy format configuration file")
            return True
            
        except Exception as e:
            self.log_warn(f"Failed to parse legacy format configuration file: {e}")
            return False
    
    def parse_config(self, config_file: str):
        """Parse configuration file"""
        if not os.path.exists(config_file):
            self.log_warn(f"Configuration file {config_file} does not exist, using default parameters")
            return
        
        config = configparser.ConfigParser()
        
        try:
            # First attempt standard INI format parsing
            config.read(config_file, encoding='utf-8')
            
            if 'analysis' in config:
                for key in self.default_config:
                    if config.has_option('analysis', key):
                        self.config[key] = config['analysis'][key]
                self.log_info(f"Successfully read standard INI configuration file: {config_file}")
            else:
                # If no analysis section, try legacy format parsing
                if not self.parse_legacy_config(config_file):
                    self.log_warn("Configuration file format is incorrect, using default parameters")
        
        except configparser.MissingSectionHeaderError:
            # Handle legacy format configuration files without section headers
            if not self.parse_legacy_config(config_file):
                self.log_warn("Configuration file format is incorrect, using default parameters")
        except Exception as e:
            self.log_warn(f"Failed to read configuration file: {e}, using default parameters")


    def check_dependencies(self):
        """Check dependencies"""
        dependencies = ['blastn', 'awk', 'sort']
        missing_deps = []
        
        for dep in dependencies:
            if not shutil.which(dep):
                missing_deps.append(dep)
        
        if missing_deps:
            self.log_error(f"Missing dependencies: {', '.join(missing_deps)}")
        
        # Check Python modules
        try:
            import pandas
            import matplotlib
            self.log_info("Python pandas and matplotlib modules are installed")
        except ImportError:
            self.log_warn("Python pandas or matplotlib modules not installed, some features may be limited")
            self.log_info("Please install: pip install pandas matplotlib")
        
        self.log_info("All dependency checks passed")
    
    def extract_random_sequences(self, reference_genome: str, extract_length: int, num_regions: int) -> Tuple[str, str]:
        """Randomly extract sequences from reference genome"""
        self.log_info(f"Randomly extracting {num_regions} region(s) of {extract_length}bp from reference genome")
        
        # Get chromosome lengths
        chrom_lengths = self.get_chromosome_lengths(reference_genome)
        if not chrom_lengths:
            self.log_error("Unable to get chromosome lengths from reference genome")
        
        # Generate random regions
        regions = []
        for i in range(num_regions):
            chrom = random.choice(list(chrom_lengths.keys()))
            max_start = chrom_lengths[chrom] - extract_length
            if max_start <= 0:
                self.log_warn(f"Chromosome {chrom} length {chrom_lengths[chrom]} is shorter than extraction length {extract_length}, skipping")
                continue
            
            start = random.randint(0, max_start)
            end = start + extract_length - 1
            regions.append((chrom, start, end, i + 1))
        
        if not regions:
            self.log_error("Unable to generate valid random regions")
        
        # Extract sequences
        output_fasta = self.output_dir / "fasta_files" / f"random_queries_{extract_length}bp.fasta"
        regions_file = self.output_dir / "fasta_files" / f"random_queries_regions.txt"
        
        with open(regions_file, 'w') as regions_f, open(output_fasta, 'w') as fasta_f:
            regions_f.write("RegionID\tChromosome\tStart\tEnd\tLength\n")
            for chrom, start, end, region_id in regions:
                regions_f.write(f"{region_id}\t{chrom}\t{start}\t{end}\t{extract_length}\n")
                
                # Use samtools to extract sequence (if available)
                if shutil.which('samtools'):
                    try:
                        result = subprocess.run([
                            'samtools', 'faidx', reference_genome,
                            f"{chrom}:{start+1}-{end+1}"  # samtools uses 1-based coordinates
                        ], capture_output=True, text=True, check=True)
                        
                        # Modify sequence header
                        lines = result.stdout.strip().split('\n')
                        if lines:
                            lines[0] = f">region_{region_id}_{chrom}_{start+1}_{end+1}"
                            fasta_f.write('\n'.join(lines) + '\n')
                        continue
                    except (FileNotFoundError, subprocess.CalledProcessError) as e:
                        self.log_warn(f"Failed to extract sequence with samtools: {e}, trying fallback method")
                
                # Fallback method: use seqkit
                if shutil.which('seqkit'):
                    try:
                        result = subprocess.run([
                            'seqkit', 'subseq', '--chr', chrom,
                            '-r', f"{start+1}:{end+1}", reference_genome
                        ], capture_output=True, text=True, check=True)
                        
                        lines = result.stdout.strip().split('\n')
                        if lines and lines[0].startswith('>'):
                            lines[0] = f">region_{region_id}_{chrom}_{start+1}_{end+1}"
                            fasta_f.write('\n'.join(lines) + '\n')
                        continue
                    except subprocess.CalledProcessError as e:
                        self.log_warn(f"Failed to extract sequence with seqkit: {e}")
                
                # Finally use Python method
                self.extract_sequences_python(reference_genome, chrom, start, end, region_id, fasta_f)
        
        self.log_info(f"Random sequences saved to: {output_fasta}")
        self.log_info(f"Region information saved to: {regions_file}")
        
        return str(output_fasta), str(regions_file)
    
    def extract_sequences_python(self, reference_genome: str, chrom: str, start: int, end: int,
                               region_id: int, fasta_file):
        """Extract sequences using Python"""
        try:
            with open(reference_genome, 'r') as f:
                in_target = False
                sequence = []
                current_chrom = ""
                
                for line in f:
                    if line.startswith('>'):
                        if in_target and sequence:
                            # Process current chromosome's sequence
                            full_sequence = ''.join(sequence)
                            if len(full_sequence) >= end:
                                extracted_seq = full_sequence[start:end+1]
                                fasta_file.write(f">region_{region_id}_{chrom}_{start+1}_{end+1}\n")
                                # 80 characters per line
                                for i in range(0, len(extracted_seq), 80):
                                    fasta_file.write(extracted_seq[i:i+80] + '\n')
                            return
                        
                        current_chrom = line[1:].split()[0].strip()
                        if current_chrom == chrom:
                            in_target = True
                            sequence = []
                        else:
                            in_target = False
                    elif in_target:
                        sequence.append(line.strip())
                
                # Process the last chromosome
                if in_target and sequence:
                    full_sequence = ''.join(sequence)
                    if len(full_sequence) >= end:
                        extracted_seq = full_sequence[start:end+1]
                        fasta_file.write(f">region_{region_id}_{chrom}_{start+1}_{end+1}\n")
                        for i in range(0, len(extracted_seq), 80):
                            fasta_file.write(extracted_seq[i:i+80] + '\n')
        
        except Exception as e:
            self.log_warn(f"Failed to extract sequence using Python method: {e}")
    
    def get_chromosome_lengths(self, fasta_file: str) -> Dict[str, int]:
        """Get chromosome lengths"""
        chrom_lengths = {}
        try:
            with open(fasta_file, 'r') as f:
                current_chrom = None
                current_length = 0
                
                for line in f:
                    if line.startswith('>'):
                        if current_chrom and current_length > 0:
                            chrom_lengths[current_chrom] = current_length
                        
                        current_chrom = line[1:].split()[0].strip()
                        current_length = 0
                    else:
                        current_length += len(line.strip())
                
                if current_chrom and current_length > 0:
                    chrom_lengths[current_chrom] = current_length
            
            self.log_info(f"Successfully read {len(chrom_lengths)} chromosome lengths")
            return chrom_lengths
        
        except Exception as e:
            self.log_warn(f"Failed to get chromosome lengths: {e}")
            return {}
    
    def read_chromosome_lengths_file(self, chr_lengths_file: str) -> Dict[str, int]:
        """Read chromosome length file (handle headers)"""
        chrom_lengths = {}
        try:
            with open(chr_lengths_file, 'r') as f:
                for i, line in enumerate(f):
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    
                    parts = line.split('\t')
                    if len(parts) >= 2:
                        chrom = parts[0].strip()
                        # Skip header row
                        if i == 0 and (chrom.lower() == 'chromosome' or chrom.lower() == 'chr' or
                                      parts[1].lower() == 'length'):
                            continue
                        try:
                            length = int(parts[1])
                            chrom_lengths[chrom] = length
                        except ValueError:
                            self.log_warn(f"Unable to parse chromosome length: {line}")
            
            self.log_info(f"Read {len(chrom_lengths)} chromosome lengths from file")
            return chrom_lengths
        except Exception as e:
            self.log_warn(f"Failed to read chromosome length file: {e}")
            return {}
    
    def generate_chromosome_lengths_file(self):
        """Generate chromosome length file"""
        reference_genome = self.config.get('reference_genome')
        if not reference_genome:
            self.log_warn("Reference genome file not provided, cannot generate chromosome length file")
            return None
        
        if not os.path.exists(reference_genome):
            self.log_warn(f"Reference genome file does not exist: {reference_genome}")
            return None
        
        # Generate chromosome length file path
        chr_lengths_file = self.output_dir / "config_files" / "chromosome_lengths.txt"
        
        # Get chromosome lengths
        chrom_lengths = self.get_chromosome_lengths(reference_genome)
        if not chrom_lengths:
            self.log_warn("Unable to get chromosome lengths from reference genome")
            return None
        
        # Save to file
        try:
            with open(chr_lengths_file, 'w') as f:
                f.write("Chromosome\tLength\n")
                for chrom, length in chrom_lengths.items():
                    f.write(f"{chrom}\t{length}\n")
            
            self.log_info(f"Chromosome length file generated: {chr_lengths_file}")
            self.config['chr_lengths'] = str(chr_lengths_file)
            self.chrom_lengths = chrom_lengths
            return str(chr_lengths_file)
        except Exception as e:
            self.log_warn(f"Failed to generate chromosome length file: {e}")
            return None
            
    def validate_parameters(self):
        """Validate parameters"""
        if not self.config['query'] and not self.config.get('reference_genome'):
            self.log_error("Must specify either query sequence file or reference genome file")
        
        if not self.config['genome_db']:
            self.log_error("Must specify BLAST database")
        
        # Validate file existence
        if self.config['query'] and not os.path.exists(self.config['query']):
            self.log_error(f"Query sequence file does not exist: {self.config['query']}")
        
        # Create output directory
        self.create_output_directory()
    
    def check_blast_database(self, db_name: str) -> bool:
        """Check BLAST database"""
        self.log_info(f"Checking BLAST database: {db_name}")
        
        # Check database files
        db_extensions = ['.nsq', '.ndb', '.nhr']
        db_exists = any(os.path.exists(db_name + ext) for ext in db_extensions)
        
        if not db_exists:
            # Check if it's a FASTA file
            if os.path.exists(db_name) and (db_name.endswith(('.fasta', '.fa', '.fna')) or
                                          os.path.exists(db_name + '.fasta') or
                                          os.path.exists(db_name + '.fa')):
                self.log_info(f"FASTA file found, attempting to create BLAST database")
                return self.create_blast_database(db_name)
            else:
                self.log_error(f"BLAST database does not exist: {db_name}")
        
        # Check database integrity
        try:
            subprocess.run(['blastdbcmd', '-db', db_name, '-info'],
                         capture_output=True, check=True)
            self.log_info("BLAST database check passed")
            return True
        except subprocess.CalledProcessError:
            self.log_warn("BLAST database may be corrupted, attempting to recreate")
            # Try to find corresponding FASTA file
            fasta_candidates = [
                db_name + '.fasta', db_name + '.fa', db_name,
                db_name.replace('.nsq', '.fasta').replace('.ndb', '.fasta').replace('.nhr', '.fasta')
            ]
            for fasta_file in fasta_candidates:
                if os.path.exists(fasta_file):
                    return self.create_blast_database(fasta_file)
            self.log_error("Cannot find corresponding FASTA file to recreate database")
    
    def create_blast_database(self, fasta_file: str) -> bool:
        """Create BLAST database"""
        if not os.path.exists(fasta_file):
            self.log_error(f"FASTA file does not exist: {fasta_file}")
        
        db_name = Path(fasta_file).stem
        
        try:
            subprocess.run([
                'makeblastdb', '-in', fasta_file, '-dbtype', 'nucl',
                '-out', db_name, '-parse_seqids'
            ], check=True)
            self.log_info(f"BLAST database created successfully: {db_name}")
            self.config['genome_db'] = db_name
            return True
        except subprocess.CalledProcessError as e:
            self.log_error(f"Failed to create BLAST database: {e}")
    
    def check_existing_results(self, args):
        """Check existing analysis result files"""
        self.log_info("Checking existing analysis results...")
        
        # Define key result files (using subdirectory paths)
        result_files = {
            'blast_output': self.output_dir / "blast_results" / f"{Path(self.config['output_prefix']).name}_blast.txt",
            'blast_sorted': self.output_dir / "blast_results" / f"{Path(self.config['output_prefix']).name}_blast_sorted.txt",
            'analysis_output': self.output_dir / "repeat_analysis" / f"{Path(self.config['output_prefix']).name}_repeat_units.txt",
            'top_repeats': self.output_dir / "repeat_analysis" / f"{Path(self.config['output_prefix']).name}_top_repeats.txt",
            'high_quality': self.output_dir / "repeat_analysis" / f"{Path(self.config['output_prefix']).name}_high_quality_repeats.txt",
            'details_output': self.output_dir / "repeat_analysis" / f"{Path(self.config['output_prefix']).name}_top_repeats_details.txt",
            'statistics_report': self.output_dir / "repeat_analysis" / f"{Path(self.config['output_prefix']).name}_statistics_report.txt"
        }
        
        # Check each file for existence
        for file_key, file_path in result_files.items():
            if file_path.exists():
                self.existing_results[file_key] = str(file_path)
                self.log_info(f"Found existing file: {file_path}")
        
        if args.force:
            self.log_info("Force re-running all steps, ignoring existing results")
            self.existing_results.clear()
        
        return len(self.existing_results) > 0
    
    def run_blast(self) -> str:
        """Run BLAST analysis"""
        # Check for existing results
        if 'blast_output' in self.existing_results:
            self.log_info(f"Using existing BLAST results: {self.existing_results['blast_output']}")
            return self.existing_results['blast_output']
        
        self.log_info("Starting BLAST analysis...")
        
        # If query sequence is empty, randomly extract from reference genome
        if not self.config['query'] and self.config.get('reference_genome'):
            if not os.path.exists(self.config['reference_genome']):
                self.log_error(f"Reference genome file does not exist: {self.config['reference_genome']}")
            
            extract_length = int(self.config['extract_length'])
            num_regions = int(self.config['extract_regions'])
            self.config['query'], regions_file = self.extract_random_sequences(
                self.config['reference_genome'], extract_length, num_regions
            )
        
        # Directly save to blast_results subdirectory
        blast_output = self.output_dir / "blast_results" / f"{Path(self.config['output_prefix']).name}_blast.txt"
        
        # Check database
        if not self.check_blast_database(self.config['genome_db']):
            self.log_error("BLAST database check failed")
        
        try:
            subprocess.run([
                'blastn', '-query', self.config['query'],
                '-db', self.config['genome_db'],
                '-out', str(blast_output),
                '-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore',
                '-evalue', self.config['e_value'],
                '-num_threads', self.config['threads'],
                '-task', 'blastn'
            ], check=True)
            
            self.log_info(f"BLAST analysis completed: {blast_output}")
            return str(blast_output)
        
        except subprocess.CalledProcessError as e:
            self.log_error(f"BLAST analysis failed: {e}")
    
    def preprocess_blast_results(self, blast_output: str) -> str:
        """Preprocess BLAST results"""
        # Check for existing results
        if 'blast_sorted' in self.existing_results:
            self.log_info(f"Using existing sorted BLAST results: {self.existing_results['blast_sorted']}")
            return self.existing_results['blast_sorted']
        
        # Directly save to blast_results subdirectory
        processed_output = self.output_dir / "blast_results" / f"{Path(self.config['output_prefix']).name}_blast_sorted.txt"
        
        self.log_info("Preprocessing BLAST results: sorting...")
        
        try:
            # Use system sort command for sorting (more efficient for large files)
            with open(processed_output, 'w') as out_file:
                subprocess.run([
                    'sort', '-k1,1', '-k7,7n', '-k8,8n', blast_output
                ], stdout=out_file, check=True)
            
            self.log_info(f"Preprocessing completed: {processed_output}")
            return str(processed_output)
        
        except subprocess.CalledProcessError as e:
            self.log_error(f"Failed to sort BLAST results: {e}")

    def analyze_repeat_units(self, blast_sorted: str) -> str:
        """Analyze repeat sequence units"""
        # Check for existing results
        if 'analysis_output' in self.existing_results:
            self.log_info(f"Using existing repeat analysis results: {self.existing_results['analysis_output']}")
            return self.existing_results['analysis_output']
        
        # Directly save to repeat_analysis subdirectory
        analysis_output = self.output_dir / "repeat_analysis" / f"{Path(self.config['output_prefix']).name}_repeat_units.txt"
        top_repeats_output = self.output_dir / "repeat_analysis" / f"{Path(self.config['output_prefix']).name}_top_repeats.txt"
        high_quality_output = self.output_dir / "repeat_analysis" / f"{Path(self.config['output_prefix']).name}_high_quality_repeats.txt"
        
        min_id = float(self.config['min_identity'])
        min_len = int(self.config['min_length'])
        min_copies = int(self.config['min_copies'])
        min_avg_id = float(self.config['min_avg_identity'])
        min_avg_len = int(self.config['min_avg_length'])
        top_n = int(self.config['top_n'])
        
        self.log_info(f"Analyzing repeat sequence units...")
        self.log_info(f"Filtering criteria: min identity {min_id}%, min length {min_len}bp, min copies {min_copies}")
        self.log_info(f"High-quality standards: avg identity ≥{min_avg_id}%, avg length ≥{min_avg_len}bp")
        
        try:
            # Read BLAST results
            df = pd.read_csv(blast_sorted, sep='\t', header=None,
                            names=['qseqid', 'sseqid', 'pident', 'length', 'mismatch',
                                  'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'])
            
            # Filter basic conditions
            df_filtered = df[(df['pident'] >= min_id) & (df['length'] >= min_len)].copy()
            
            if df_filtered.empty:
                self.log_warn("No BLAST results meeting basic conditions")
                # Create empty output file
                with open(analysis_output, 'w') as f:
                    f.write("RepeatUnitID\tQueryID\tQueryStart\tQueryEnd\tCopyNumber\tAvgIdentity\tAvgLength\tTotalLength\tAvgEvalue\tAvgScore\tQuality\n")
                return str(analysis_output)
            
            # Define repeat units
            df_filtered.loc[:, 'unit_id'] = df_filtered['qseqid'] + ':' + df_filtered['qstart'].astype(str) + ':' + df_filtered['qend'].astype(str)
            
            # Aggregate repeat unit statistics
            repeat_stats = df_filtered.groupby('unit_id').agg({
                'pident': ['count', 'mean'],
                'length': 'mean',
                'evalue': 'mean',
                'bitscore': 'mean',
                'qseqid': 'first',
                'qstart': 'first',
                'qend': 'first'
            }).round(3)
            
            # Rename columns
            repeat_stats.columns = ['copy_number', 'avg_identity', 'avg_length', 'avg_evalue', 'avg_score', 'query_id', 'query_start', 'query_end']
            repeat_stats['total_length'] = repeat_stats['copy_number'] * repeat_stats['avg_length']
            
            # Add quality label - first determine if high quality
            repeat_stats['quality'] = 'Standard'
            high_quality_mask = (repeat_stats['avg_identity'] >= min_avg_id) & (repeat_stats['avg_length'] >= min_avg_len)
            repeat_stats.loc[high_quality_mask, 'quality'] = 'HighQuality'
            
            # Filter by minimum copy number
            repeat_stats = repeat_stats[repeat_stats['copy_number'] >= min_copies]
            
            # Reset index to save unit_id
            repeat_stats_reset = repeat_stats.reset_index()
            
            # Reorder columns
            column_order = ['unit_id', 'query_id', 'query_start', 'query_end', 'copy_number',
                          'avg_identity', 'avg_length', 'total_length', 'avg_evalue', 'avg_score', 'quality']
            repeat_stats_reset = repeat_stats_reset[column_order]
            
            # Get high-quality repeat sequences (for top selection and separate saving)
            high_quality_df = repeat_stats_reset[repeat_stats_reset['quality'] == 'HighQuality'].copy()
            
            # Statistics
            total_units = len(repeat_stats)
            high_quality_units = len(high_quality_df)
            
            self.log_info(f"Found {total_units} repeat units meeting conditions (copy number ≥{min_copies})")
            self.log_info(f"Of these, {high_quality_units} are high-quality repeat units")
            
            if total_units == 0:
                self.log_warn("No repeat units meeting conditions found")
                # Save debug information
                debug_file = self.output_dir / "repeat_analysis" / f"{Path(self.config['output_prefix']).name}_debug.txt"
                debug_stats = df_filtered.groupby('unit_id').agg({
                    'pident': ['count', 'mean']
                }).round(3)
                debug_stats.columns = ['copy_number', 'avg_identity']
                debug_stats = debug_stats.sort_values('copy_number', ascending=False).head(20)
                debug_stats.to_csv(debug_file, sep='\t')
                self.log_info(f"Debug information saved to: {debug_file}")
            
            # Save all repeat units
            repeat_stats_reset.to_csv(analysis_output, sep='\t', index=False)
            
            # Save high-quality repeat sequences
            if high_quality_units > 0:
                high_quality_df.to_csv(high_quality_output, sep='\t', index=False)
                self.log_info(f"High-quality repeat sequences saved to: {high_quality_output}")
                
                # Select top N from high-quality repeat sequences (based on copy number)
                if high_quality_units > top_n:
                    top_n_df = high_quality_df.sort_values('copy_number', ascending=False).head(top_n)
                    self.log_info(f"Selecting top {top_n} from high-quality repeat sequences (based on copy number)")
                else:
                    top_n_df = high_quality_df.sort_values('copy_number', ascending=False)
                    self.log_info(f"Number of high-quality repeat sequences ({high_quality_units}) is less than top_n ({top_n}), saving all high-quality repeats")
                
                top_n_df.to_csv(top_repeats_output, sep='\t', index=False)
            else:
                self.log_warn("No high-quality repeat sequences, skipping top repeats file generation")
                # Save an empty top file
                with open(top_repeats_output, 'w') as f:
                    f.write("RepeatUnitID\tQueryID\tQueryStart\tQueryEnd\tCopyNumber\tAvgIdentity\tAvgLength\tTotalLength\tAvgEvalue\tAvgScore\tQuality\n")
            
            self.log_info(f"Repeat unit analysis completed: {analysis_output}")
            return str(analysis_output)
        
        except Exception as e:
            self.log_error(f"Failed to analyze repeat units: {e}")
    
    def extract_top_repeats_details(self, blast_sorted: str, analysis_output: str) -> str:
        """Extract detailed information for top repeat sequences (now only from high-quality repeats)"""
        # Check for existing results
        if 'details_output' in self.existing_results:
            self.log_info(f"Using existing detailed information: {self.existing_results['details_output']}")
            return self.existing_results['details_output']
        
        # Directly save to repeat_analysis subdirectory
        details_output = self.output_dir / "repeat_analysis" / f"{Path(self.config['output_prefix']).name}_top_repeats_details.txt"
        top_n = int(self.config['top_n'])
        
        self.log_info(f"Extracting detailed information for top {top_n} high-quality repeat sequences...")
        
        try:
            # Read high-quality repeat sequences file
            high_quality_file = self.output_dir / "repeat_analysis" / f"{Path(self.config['output_prefix']).name}_high_quality_repeats.txt"
            
            if not high_quality_file.exists():
                self.log_warn("High-quality repeat sequences file does not exist, skipping detailed information extraction")
                # Create empty details file
                with open(details_output, 'w') as f:
                    f.write("unit_id\tqseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstatus\n")
                return str(details_output)
            
            high_quality_df = pd.read_csv(high_quality_file, sep='\t')
            
            if high_quality_df.empty:
                self.log_warn("High-quality repeat sequences are empty, skipping detailed information extraction")
                # Create empty details file
                with open(details_output, 'w') as f:
                    f.write("unit_id\tqseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstatus\n")
                return str(details_output)
            
            # Select top N from high-quality repeat sequences (based on copy number)
            top_units_df = high_quality_df.sort_values('copy_number', ascending=False).head(top_n)
            top_units = top_units_df['unit_id'].tolist()
            
            self.log_info(f"Selected {len(top_units)} top repeat units from high-quality sequences")
            
            # Read BLAST results
            df = pd.read_csv(blast_sorted, sep='\t', header=None,
                            names=['qseqid', 'sseqid', 'pident', 'length', 'mismatch',
                                  'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'])
            
            # Add repeat unit ID
            df['unit_id'] = df['qseqid'] + ':' + df['qstart'].astype(str) + ':' + df['qend'].astype(str)
            
            # Filter top repeat units
            df_top = df[df['unit_id'].isin(top_units)].copy()
            
            # Add filter status
            min_id = float(self.config['min_identity'])
            min_len = int(self.config['min_length'])
            df_top['status'] = 'PASS'
            df_top.loc[(df_top['pident'] < min_id) | (df_top['length'] < min_len), 'status'] = 'FILTERED'
            
            # Reorder columns
            column_order = ['unit_id', 'qseqid', 'sseqid', 'pident', 'length', 'mismatch',
                          'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'status']
            df_top = df_top[column_order]
            
            # Save detailed information
            df_top.to_csv(details_output, sep='\t', index=False)
            
            # Statistics
            total_details = len(df_top)
            passed_details = len(df_top[df_top['status'] == 'PASS'])
            filtered_details = len(df_top[df_top['status'] == 'FILTERED'])
            
            self.log_info(f"Detailed information statistics: total entries {total_details}, passed {passed_details}, filtered {filtered_details}")
            self.log_info(f"Detailed information extraction completed: {details_output}")
            
            return str(details_output)
        
        except Exception as e:
            self.log_error(f"Failed to extract detailed information: {e}")

    def generate_high_quality_fasta(self, analysis_output: str) -> str:
        """Generate FASTA file for high-quality repeat sequences"""
        if not self.config['query']:
            self.log_info("Using random extracted sequences, skipping high-quality repeat FASTA generation")
            return ""
        
        # Directly save to fasta_files subdirectory
        fasta_output = self.output_dir / "fasta_files" / f"{Path(self.config['output_prefix']).name}_high_quality_repeats.fasta"
        
        # Check for existing results
        if fasta_output.exists():
            self.log_info(f"Using existing FASTA file: {fasta_output}")
            return str(fasta_output)
        
        self.log_info("Generating FASTA file for high-quality repeat sequences...")
        
        try:
            # Read high-quality repeat sequences
            repeat_stats = pd.read_csv(analysis_output, sep='\t')
            high_quality = repeat_stats[repeat_stats['quality'] == 'HighQuality']
            
            if high_quality.empty:
                self.log_warn("No high-quality repeat sequences, skipping FASTA file generation")
                return ""
            
            # Extract sequence region information
            regions = []
            for _, row in high_quality.iterrows():
                unit_id = row['unit_id']
                query_id = row['query_id']
                start = int(row['query_start'])
                end = int(row['query_end'])
                regions.append((query_id, start, end, unit_id))
            
            # Use seqkit to extract sequences
            if shutil.which('seqkit'):
                self.extract_with_seqkit(regions, str(fasta_output))
            else:
                self.extract_with_python(regions, str(fasta_output))
            
            self.log_info(f"High-quality repeat FASTA file generation completed: {fasta_output}")
            return str(fasta_output)
        
        except Exception as e:
            self.log_error(f"Failed to generate FASTA file: {e}")
    
    def extract_with_seqkit(self, regions: List[Tuple], fasta_output: str):
        """Extract sequences using seqkit"""
        with open(fasta_output, 'w') as out_f:
            for query_id, start, end, unit_id in regions:
                try:
                    # Use seqkit to extract sequence
                    result = subprocess.run([
                        'seqkit', 'subseq', '--chr', query_id,
                        '-r', f"{start}:{end}", self.config['query']
                    ], capture_output=True, text=True, check=True)
                    
                    # Modify sequence header
                    lines = result.stdout.strip().split('\n')
                    if lines and lines[0].startswith('>'):
                        lines[0] = f">{unit_id}_{start}-{end}"
                        out_f.write('\n'.join(lines) + '\n')
                
                except subprocess.CalledProcessError as e:
                    self.log_warn(f"Failed to extract sequence with seqkit: {unit_id}, trying fallback method")
                    self.extract_single_sequence_python(query_id, start, end, unit_id, out_f)
    
    def extract_with_python(self, regions: List[Tuple], fasta_output: str):
        """Extract sequences using Python"""
        with open(fasta_output, 'w') as out_f:
            for query_id, start, end, unit_id in regions:
                self.extract_single_sequence_python(query_id, start, end, unit_id, out_f)
    
    def extract_single_sequence_python(self, query_id: str, start: int, end: int, unit_id: str, out_file):
        """Extract a single sequence using Python method"""
        try:
            with open(self.config['query'], 'r') as f:
                in_target = False
                sequence = []
                
                for line in f:
                    if line.startswith('>'):
                        current_id = line[1:].split()[0].strip()
                        if current_id == query_id:
                            in_target = True
                            sequence = []
                        else:
                            in_target = False
                    elif in_target:
                        sequence.append(line.strip())
                
                if sequence:
                    full_sequence = ''.join(sequence)
                    if len(full_sequence) >= end:
                        extracted_seq = full_sequence[start-1:end]  # Convert to 0-based
                        out_file.write(f">{unit_id}_{start}-{end}\n")
                        # 80 characters per line
                        for i in range(0, len(extracted_seq), 80):
                            out_file.write(extracted_seq[i:i+80] + '\n')
                    else:
                        self.log_warn(f"Sequence length insufficient: {query_id}, requires positions {start}-{end}, but sequence length is only {len(full_sequence)}")
                else:
                    self.log_warn(f"Sequence not found: {query_id}")
        
        except Exception as e:
            self.log_warn(f"Failed to extract sequence {unit_id}: {e}")
    
    def generate_statistics_report(self, details_output: str):
        """Generate statistics report"""
        # Check for existing results
        if 'statistics_report' in self.existing_results:
            self.log_info(f"Using existing statistics report: {self.existing_results['statistics_report']}")
            return
        
        # Directly save to repeat_analysis subdirectory
        report_output = self.output_dir / "repeat_analysis" / f"{Path(self.config['output_prefix']).name}_statistics_report.txt"
        
        self.log_info("Generating statistics report...")
        
        try:
            df = pd.read_csv(details_output, sep='\t')
            
            # Generate statistics report
            with open(report_output, 'w') as f:
                f.write("=== Repeat Sequence Analysis Statistics Report ===\n")
                f.write(f"Filtering criteria:\n")
                f.write(f"  Minimum identity: {self.config['min_identity']}%\n")
                f.write(f"  Minimum length: {self.config['min_length']} bp\n")
                f.write(f"  Minimum copy number: {self.config['min_copies']}\n")
                f.write(f"  High-quality standard: avg identity ≥{self.config['min_avg_identity']}%, avg length ≥{self.config['min_avg_length']}bp\n\n")
                
                f.write("Repeat Unit Statistics:\n")
                
                # Aggregate by repeat unit
                unit_stats = df.groupby('unit_id').agg({
                    'pident': ['count', 'mean'],
                    'length': 'mean'
                }).round(2)
                
                unit_stats.columns = ['total_copies', 'avg_identity', 'avg_length']
                unit_stats['passed_copies'] = df[df['status'] == 'PASS'].groupby('unit_id').size()
                unit_stats['filtered_copies'] = unit_stats['total_copies'] - unit_stats['passed_copies']
                unit_stats = unit_stats.fillna(0)
                
                # Write statistics results
                f.write("RepeatUnitID\tTotalCopies\tPassedCopies\tFilteredCopies\tAvgIdentity\tAvgLength\n")
                for unit_id, row in unit_stats.iterrows():
                    f.write(f"{unit_id}\t{int(row['total_copies'])}\t{int(row['passed_copies'])}\t{int(row['filtered_copies'])}\t{row['avg_identity']}\t{row['avg_length']:.1f}\n")
            
            self.log_info(f"Statistics report generation completed: {report_output}")
        
        except Exception as e:
            self.log_warn(f"Failed to generate statistics report: {e}")

    def load_chromosome_lengths(self):
        """Load chromosome length data"""
        # If chromosome length file is provided in configuration
        if self.config['chr_lengths'] and os.path.exists(self.config['chr_lengths']):
            self.chrom_lengths = self.read_chromosome_lengths_file(self.config['chr_lengths'])
            if self.chrom_lengths:
                self.log_info(f"Loaded {len(self.chrom_lengths)} chromosome lengths from configuration file")
                # Debug: show first few chromosome IDs
                self.log_info(f"Chromosome ID examples: {list(self.chrom_lengths.keys())[:10]}")
                return True
            else:
                self.log_warn("Chromosome length file in configuration is invalid")
        
        # If not provided or invalid, try to generate from reference genome
        if self.config.get('reference_genome') and os.path.exists(self.config['reference_genome']):
            self.log_info("Attempting to generate chromosome length file from reference genome...")
            chr_lengths_file = self.generate_chromosome_lengths_file()
            if chr_lengths_file and os.path.exists(chr_lengths_file):
                self.chrom_lengths = self.read_chromosome_lengths_file(chr_lengths_file)
                if self.chrom_lengths:
                    self.log_info(f"Generated {len(self.chrom_lengths)} chromosome lengths from reference genome")
                    self.log_info(f"Chromosome ID examples: {list(self.chrom_lengths.keys())[:10]}")
                    return True
        
        self.log_warn("Unable to obtain chromosome length data, skipping distribution plot generation")
        return False
    
    def generate_distribution_plots(self, details_output: str, analysis_output: str):
        """Generate chromosome distribution plots (only for top high-quality repeat sequences)"""
        # Check if visualization results already exist
        vis_dir = self.output_dir / "visualization"
        
        # Load chromosome length data
        if not self.load_chromosome_lengths():
            return
        
        self.log_info("Generating chromosome distribution plots for top high-quality repeat units...")
        
        try:
            # Read top repeat sequences file
            top_repeats_file = self.output_dir / "repeat_analysis" / f"{Path(self.config['output_prefix']).name}_top_repeats.txt"
            
            if not top_repeats_file.exists():
                self.log_warn("Top repeat sequences file does not exist, skipping distribution plot generation")
                return
            
            top_repeats_df = pd.read_csv(top_repeats_file, sep='\t')
            
            if top_repeats_df.empty:
                self.log_warn("Top repeat sequences file is empty, skipping distribution plot generation")
                return
            
            top_units = top_repeats_df['unit_id'].tolist()
            
            # Read detailed information
            df_details = pd.read_csv(details_output, sep='\t')
            
            self.log_info(f"Generating distribution plots for {len(top_units)} high-quality top repeat sequences")
            
            # Generate distribution plots for all top high-quality repeat sequences
            plotted_count = 0
            for unit_id in top_units:
                try:
                    self.generate_single_distribution_plot(unit_id, df_details)
                    plotted_count += 1
                except Exception as e:
                    self.log_warn(f"Failed to generate distribution plot for repeat unit {unit_id}: {e}")
                    continue
            
            self.log_info(f"Successfully generated distribution plots for {plotted_count} repeat sequences")
            
            # Generate overall distribution plot (only for top high-quality repeat sequences)
            self.generate_overall_distribution_plot(df_details, top_units)
            
        except Exception as e:
            self.log_warn(f"Failed to generate distribution plots: {e}")
    
    def generate_single_distribution_plot(self, unit_id: str, df_details: pd.DataFrame):
        """Generate distribution plot for a single repeat unit"""
        try:
            # Filter data for current repeat unit
            unit_data = df_details[df_details['unit_id'] == unit_id]
            unit_data = unit_data[unit_data['status'] == 'PASS']
            
            if unit_data.empty:
                self.log_warn(f"Repeat unit {unit_id} has no data that passed filtering, skipping plot generation")
                return
            
            # Clean chromosome IDs: remove 'gb|' prefix and '|' suffix or 'ref|' prefix and '|' suffix
            def clean_chrom_id(chrom_id):
                # Remove common FASTA header formats
                chrom_str = str(chrom_id)
                if chrom_str.startswith('ref|'):
                    chrom_str = chrom_str[4:]
                if '|' in chrom_str:
                    chrom_str = chrom_str.split('|')[0]
                return chrom_str
            
            # Create cleaned chromosome ID column
            unit_data = unit_data.copy()
            unit_data['clean_sseqid'] = unit_data['sseqid'].apply(clean_chrom_id)
            
            # Check chromosome ID matching
            unmatched_chroms = set(unit_data['clean_sseqid']) - set(self.chrom_lengths.keys())
            matched_chroms = set(unit_data['clean_sseqid']) & set(self.chrom_lengths.keys())
            
            if unmatched_chroms:
                self.log_warn(f"Repeat unit {unit_id}: {len(unmatched_chroms)} chromosome IDs did not match")
            
            # Safely handle unit_id (replace special characters)
            safe_unit_id = unit_id.replace(':', '_').replace('/', '_').replace('\\', '_')
            
            # Directly save to visualization subdirectory
            plot_prefix = self.output_dir / "visualization" / f"{Path(self.config['output_prefix']).name}_{safe_unit_id}_distribution"
            
            # Create figure
            fig, ax = plt.subplots(figsize=(16, 10))
            
            # Sort chromosomes (based on chromosome length file)
            chromosomes = sorted(self.chrom_lengths.keys(),
                               key=lambda x: int(x.replace('chr', '').replace('Chr', '').replace('CHR', ''))
                               if x.replace('chr', '').replace('Chr', '').replace('CHR', '').isdigit() else x)
            
            # Count hits per chromosome
            chrom_hits = {}
            plotted_hits = 0
            for chrom in chromosomes:
                # Match cleaned chromosome ID
                chrom_data = unit_data[unit_data['clean_sseqid'] == chrom]
                chrom_hits[chrom] = len(chrom_data)
            
            # Plot each chromosome
            for i, chrom in enumerate(chromosomes):
                if chrom not in self.chrom_lengths:
                    continue
                    
                chrom_length = self.chrom_lengths[chrom]
                chrom_length_mb = chrom_length / 1000000  # Convert to Mb
                
                # Draw chromosome background
                ax.hlines(y=i, xmin=0, xmax=chrom_length_mb,
                          color='lightgray', linewidth=6, alpha=0.7)
                
                # Draw repeat sequence positions on this chromosome
                chrom_data = unit_data[unit_data['clean_sseqid'] == chrom]
                
                for _, hit in chrom_data.iterrows():
                    # Check if coordinates are within chromosome range
                    start_pos = min(hit['sstart'], hit['send'])
                    end_pos = max(hit['sstart'], hit['send'])
                    
                    if start_pos < 0 or end_pos > chrom_length:
                        continue
                    
                    identity = hit['pident']
                    if identity >= 90:
                        color = 'red'
                        alpha = 0.9
                        linewidth = 2
                    elif identity >= 80:
                        color = 'green'
                        alpha = 0.8
                        linewidth = 1.5
                    else:
                        color = 'blue'
                        alpha = 0.7
                        linewidth = 1
                    
                    # Use start position (convert to Mb)
                    pos = start_pos / 1000000
                    
                    # Draw vertical line segment
                    ax.vlines(x=pos, ymin=i-0.35, ymax=i+0.35,
                              color=color, linewidth=linewidth, alpha=alpha)
                    plotted_hits += 1
            
            # Beautify plot
            ax.set_yticks(range(len(chromosomes)))
            ax.set_yticklabels(chromosomes)
            ax.set_xlabel('Genomic Position (Mb)', fontsize=12)
            ax.set_ylabel('Chromosome', fontsize=12)
            ax.set_title(f'Distribution of Repeat Unit: {unit_id}\nTotal Hits: {len(unit_data)}, Plotted: {plotted_hits}',
                        fontsize=14, fontweight='bold')
            
            # Calculate maximum chromosome length
            max_chrom_length = max(self.chrom_lengths.values()) if self.chrom_lengths else 1
            max_chrom_length_mb = max_chrom_length / 1000000
            
            # Set axis limits
            ax.set_xlim(0, max_chrom_length_mb * 1.05)
            ax.set_ylim(-0.5, len(chromosomes) - 0.5)
            
            # Add grid
            ax.grid(True, axis='x', alpha=0.3, linestyle='--')
            
            # Add legend
            legend_elements = [
                Patch(facecolor='red', label='Identity ≥90%'),
                Patch(facecolor='green', label='Identity 80-90%'),
                Patch(facecolor='blue', label='Identity <80%')
            ]
            ax.legend(handles=legend_elements, loc='upper right')
            
            # Save plot
            plt.tight_layout()
            plt.savefig(f"{plot_prefix}.png", dpi=300, bbox_inches='tight')
            plt.savefig(f"{plot_prefix}.pdf", bbox_inches='tight')
            plt.close()
            
            self.log_info(f"Successfully generated distribution plot: {plot_prefix}.png, plotted {plotted_hits} hits")
            
        except Exception as e:
            self.log_warn(f"Failed to generate single distribution plot {unit_id}: {e}")

    def generate_overall_distribution_plot(self, df_details: pd.DataFrame, top_units: List[str]):
        """Generate overall distribution plot (only for top high-quality repeat sequences)"""
        try:
            # Filter passed data and only include top repeat sequences
            passed_data = df_details[(df_details['status'] == 'PASS') & 
                                   (df_details['unit_id'].isin(top_units))]
            
            if passed_data.empty:
                self.log_warn("No passed data for top repeat sequences, skipping overall distribution plot generation")
                return
            
            self.log_info(f"Generating overall distribution plot for high-quality top repeat sequences, total {len(passed_data)} hits")
            
            # Function to clean chromosome IDs
            def clean_chrom_id(chrom_id):
                chrom_str = str(chrom_id)
                if chrom_str.startswith('ref|'):
                    chrom_str = chrom_str[4:]
                if '|' in chrom_str:
                    chrom_str = chrom_str.split('|')[0]
                return chrom_str
            
            # Create cleaned chromosome ID column
            passed_data = passed_data.copy()
            passed_data['clean_sseqid'] = passed_data['sseqid'].apply(clean_chrom_id)
            
            # Debug: show chromosome ID matching status
            unique_chroms = passed_data['clean_sseqid'].unique()
            self.log_info(f"Number of unique cleaned chromosome IDs: {len(unique_chroms)}")
            self.log_info(f"Number of IDs in chromosome length file: {len(self.chrom_lengths)}")
            
            # Check unmatched chromosome IDs
            unmatched = set(unique_chroms) - set(self.chrom_lengths.keys())
            matched = set(unique_chroms) & set(self.chrom_lengths.keys())
            
            if unmatched:
                self.log_warn(f"{len(unmatched)} chromosome IDs not found in length file:")
                for chrom in list(unmatched)[:10]:
                    self.log_warn(f"  {chrom}")
                if len(unmatched) > 10:
                    self.log_warn(f"  ... and {len(unmatched)-10} more")
            
            self.log_info(f"Successfully matched {len(matched)} chromosome IDs")
            
            # Directly save to visualization subdirectory
            plot_prefix = self.output_dir / "visualization" / f"{Path(self.config['output_prefix']).name}_overall_distribution"
            
            # Create figure
            fig, ax = plt.subplots(figsize=(16, 10))
            
            # Sort chromosomes (based on chromosome length file)
            chromosomes = sorted(self.chrom_lengths.keys(),
                               key=lambda x: int(x.replace('chr', '').replace('Chr', '').replace('CHR', ''))
                               if x.replace('chr', '').replace('Chr', '').replace('CHR', '').isdigit() else x)
            
            # Statistics
            total_hits = 0
            plotted_hits = 0
            for chrom in chromosomes:
                chrom_data = passed_data[passed_data['clean_sseqid'] == chrom]
                hit_count = len(chrom_data)
                if hit_count > 0:
                    self.log_info(f"  Chromosome {chrom}: {hit_count} hits")
                    total_hits += hit_count
            
            # Plot each chromosome
            for i, chrom in enumerate(chromosomes):
                if chrom not in self.chrom_lengths:
                    continue
                    
                chrom_length = self.chrom_lengths[chrom]
                chrom_length_mb = chrom_length / 1000000
                
                # Draw chromosome background
                ax.hlines(y=i, xmin=0, xmax=chrom_length_mb,
                          color='lightgray', linewidth=4, alpha=0.7)
                
                # Draw repeat sequence positions on this chromosome
                chrom_data = passed_data[passed_data['clean_sseqid'] == chrom]
                
                for _, hit in chrom_data.iterrows():
                    # Check if coordinates are within chromosome range
                    start_pos = min(hit['sstart'], hit['send'])
                    end_pos = max(hit['sstart'], hit['send'])
                    
                    if start_pos < 0 or end_pos > chrom_length:
                        self.log_warn(f"Coordinates out of range: chromosome {chrom}, positions {start_pos}-{end_pos}, chromosome length {chrom_length}")
                        continue
                    
                    identity = hit['pident']
                    if identity >= 90:
                        color = 'red'
                        alpha = 0.6
                    elif identity >= 80:
                        color = 'green'
                        alpha = 0.5
                    else:
                        color = 'blue'
                        alpha = 0.4
                    
                    # Use start position (convert to Mb)
                    pos = start_pos / 1000000
                    
                    # Draw vertical line segment
                    ax.vlines(x=pos, ymin=i-0.25, ymax=i+0.25,
                              color=color, linewidth=0.8, alpha=alpha)
                    plotted_hits += 1
            
            self.log_info(f"Total {total_hits} hits, successfully plotted {plotted_hits}")
            
            if plotted_hits == 0:
                self.log_warn("No hits were successfully plotted!")
                # Add debug information to the plot
                ax.text(0.5, 0.5, 'No hits plotted!\nPossible reasons:\n1. Chromosome ID mismatch\n2. Coordinates out of range',
                       transform=ax.transAxes, fontsize=14, ha='center', va='center',
                       bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
            
            # Beautify plot
            ax.set_yticks(range(len(chromosomes)))
            ax.set_yticklabels(chromosomes)
            ax.set_xlabel('Genomic Position (Mb)', fontsize=12)
            ax.set_ylabel('Chromosome', fontsize=12)
            ax.set_title(f'Genome-wide Distribution of High-Quality Top Repeat Sequences\nTotal Hits: {total_hits}, Plotted: {plotted_hits}',
                        fontsize=14, fontweight='bold')
            
            # Calculate maximum chromosome length
            max_chrom_length = max(self.chrom_lengths.values()) if self.chrom_lengths else 1
            max_chrom_length_mb = max_chrom_length / 1000000
            
            # Set axis limits
            ax.set_xlim(0, max_chrom_length_mb * 1.05)
            ax.set_ylim(-0.5, len(chromosomes) - 0.5)
            
            # Add grid
            ax.grid(True, axis='x', alpha=0.3, linestyle='--')
            
            # Add legend
            legend_elements = [
                Patch(facecolor='red', label='Identity ≥90%'),
                Patch(facecolor='green', label='Identity 80-90%'),
                Patch(facecolor='blue', label='Identity <80%')
            ]
            ax.legend(handles=legend_elements, loc='upper right')
            
            # Save plot
            plt.tight_layout()
            plt.savefig(f"{plot_prefix}.png", dpi=300, bbox_inches='tight')
            plt.savefig(f"{plot_prefix}.pdf", bbox_inches='tight')
            plt.close()
            
            self.log_info(f"Successfully generated overall distribution plot: {plot_prefix}.png")
            
        except Exception as e:
            self.log_warn(f"Failed to generate overall distribution plot: {e}")
            import traceback
            traceback.print_exc()
    
    def organize_results(self):
        """Organize result files (now mainly handles moving log files and configuration files)"""
        self.log_info("Organizing result files...")
        
        # Move log file
        try:
            if os.path.exists('repeat_analysis.log'):
                log_target = self.output_dir / 'repeat_analysis.log'
                shutil.move('repeat_analysis.log', str(log_target))
                self.log_info(f"Moved log file: repeat_analysis.log -> {self.output_dir}")
        except Exception as e:
            self.log_warn(f"Failed to move log file: {e}")
        
        # Copy configuration file
        try:
            args = self.parse_arguments()
            config_file = Path(args.config_file)
            if config_file.exists():
                config_target = self.output_dir / "config_files" / config_file.name
                config_target.parent.mkdir(exist_ok=True, parents=True)
                shutil.copy2(str(config_file), str(config_target))
                self.log_info(f"Copied configuration file: {config_file} -> {config_target}")
        except Exception as e:
            self.log_warn(f"Failed to copy configuration file: {e}")
        
        # Check for leftover files in working directory (only show warnings, don't move)
        try:
            import glob
            leftover_files = []
            for file_pattern in ['*_blast.txt', '*_blast_sorted.txt', '*_repeat_units.txt',
                               '*_top_repeats.txt', '*_high_quality_repeats.txt', 
                               '*_top_repeats_details.txt', '*_statistics_report.txt',
                               '*_debug.txt', '*_distribution.png', '*_distribution.pdf',
                               '*_high_quality_repeats.fasta']:
                leftover_files.extend(glob.glob(file_pattern))
            
            if leftover_files:
                self.log_warn(f"Found {len(leftover_files)} files left in working directory, these files should have been saved to output directory")
                for file in leftover_files[:5]:  # Show only first 5
                    self.log_warn(f"  - {file}")
                if len(leftover_files) > 5:
                    self.log_warn(f"  ... and {len(leftover_files)-5} more files")
        except Exception as e:
            self.log_warn(f"Failed to check for leftover files: {e}")
        
        self.log_info(f"All result files have been saved to directory: {self.output_dir}")
        self.log_info(f"Subdirectory structure:")
        self.log_info(f"  - blast_results: BLAST related results")
        self.log_info(f"  - repeat_analysis: Repeat sequence analysis results")
        self.log_info(f"  - visualization: Distribution plot files")
        self.log_info(f"  - fasta_files: FASTA sequence files")
        self.log_info(f"  - config_files: Configuration files")
    
    def run(self):
        """Run the main analysis pipeline"""
        args = self.parse_arguments()
        
        if args.create_config:
            self.create_config_template()
            return
        
        self.log_info("Starting repeat sequence analysis")
        
        # Check dependencies
        self.check_dependencies()
        
        # Parse configuration
        self.parse_config(args.config_file)
        
        # Validate parameters
        self.validate_parameters()
        
        # Check existing results
        has_existing_results = self.check_existing_results(args)
        
        if has_existing_results and not args.force:
            self.log_info("Detected existing analysis results, will continue using existing results")
            self.log_info("To re-run all steps, use the --force parameter")
        else:
            self.log_info("No existing analysis results found or force re-run requested, starting full analysis pipeline")
        
        # Run BLAST
        blast_output = self.run_blast()
        
        # Preprocess BLAST results
        blast_sorted = self.preprocess_blast_results(blast_output)
        
        # Analyze repeat units
        analysis_output = self.analyze_repeat_units(blast_sorted)
        
        # Extract detailed information
        details_output = self.extract_top_repeats_details(blast_sorted, analysis_output)
        
        # Generate high-quality repeat FASTA
        fasta_output = self.generate_high_quality_fasta(analysis_output)
        
        # Generate statistics report
        self.generate_statistics_report(details_output)
        
        # Generate distribution plots
        self.generate_distribution_plots(details_output, analysis_output)
        
        # Organize result files
        self.organize_results()
        
        self.log_info("Analysis completed!")
        self.log_info(f"All result files have been saved to directory: {self.output_dir}")
        self.log_info(f"Subdirectory structure:")
        self.log_info(f"  - blast_results: BLAST related results")
        self.log_info(f"  - repeat_analysis: Repeat sequence analysis results")
        self.log_info(f"  - visualization: Distribution plot files")
        self.log_info(f"  - fasta_files: FASTA sequence files")
        self.log_info(f"  - config_files: Configuration files")
        self.log_info(f"Note: Next run will automatically use existing results. To re-run, add the --force parameter")

def main():
    """Main function"""
    analyzer = RepeatAnalyzer()
    analyzer.run()

if __name__ == '__main__':
    main()