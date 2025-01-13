# FASTA Cleaner with fixed AT ratios

## Overview

The FASTA Cleaner is a Python tool designed to analyze and filter DNA sequences based on AT content while generating consensus sequences. The tool processes aligned FASTA files and creates separate outputs for sequences that meet specified AT content criteria. This tool was developed to remove contaminant sequences (e.g. Fungi) from alignments of insect mitochondrial DNA which typically have a bias towards AT content.

## Features

The processor provides several key capabilities:
- Filters sequences based on customizable AT content thresholds
- Generates consensus sequences from filtered alignments
- Preserves alignment structure in outputs
- Processes multiple FASTA files in batch
- Outputs sequences in single-line format for computational processing

## Installation

Requirements:
- Python 3.6 or higher
- BioPython library

Install the required dependency using pip:
```bash
pip install biopython
```

## Usage

Basic command structure:
```bash
python fasta_cleaner.py -i INPUT_DIR -o OUTPUT_DIR [-t THRESHOLD] [-c CONSENSUS_THRESHOLD]
```

### Arguments

Required:
- `-i, --input_dir`: Directory containing input FASTA files
- `-o, --output_dir`: Directory where output files will be saved

Optional:
- `-t, --threshold`: AT content threshold (default: 0.5)
  - Range: 0.0-1.0
  - Sequences with AT content ≥ threshold are kept
- `-c, --consensus_threshold`: Consensus calling threshold (default: 0.7)
  - Range: 0.0-1.0
  - Minimum frequency required for a base to be included in consensus

### Output Files

For each input file named 'example.fasta', the tool generates:
1. `example_cleaned.fasta`
   - Contains sequences meeting the AT threshold
   - Full sequences on single lines
   - Maintains alignment structure
   
2. `example_removed.fasta`
   - Contains sequences below the AT threshold
   - Same format as cleaned file
   - Useful for quality control and verification
   
3. `example_cleaned_consensus.fasta`
   - Consensus sequence derived from cleaned sequences
   - Uses specified consensus threshold
   - Marks ambiguous positions with 'X'

4. `example_report.csv`
   - One row for each sequence in the original file
   - Sequence name
   - Sequence length (excluding gaps)
   - AT ratio
   - Status (kept / removed)
     
## Example Usage

Basic usage with default thresholds:
```bash
python fasta_cleaner.py -i /path/to/fasta/files -o /path/to/results
```

Custom thresholds for more stringent filtering:
```bash
python fasta_cleaner.py -i /path/to/fasta/files -o /path/to/results -t 0.6 -c 0.8
```

## Considerations

1. AT Content Threshold:
   - Higher values (e.g., 0.6) are more stringent
   - Consider your organism's typical AT content
   - Some regions of a gene may have higher GC content
   - Test different thresholds on a subset of data first

2. Consensus Threshold:
   - Higher values (e.g., 0.8) give more conservative consensus
   - Lower values (e.g., 0.5) mask more variable positions
   - Adjust based on sequence variability and analysis needs

3. Input Files:
   - Must be aligned FASTA format
   - Supports .fasta, .fas, and .fa extensions
   - Maintains gap positions in alignments
