# FASTA Cleaner with contamination detection

A robust Python tool for cleaning FASTA sequence alignments by identifying and removing potential contaminant sequences through protein-level comparisons.

## Features

- Processes multiple FASTA files in batch
- Translates DNA sequences in all three reading frames
- Performs protein-level alignments using BLOSUM62 scoring matrix
- Compares sequences against both target and contaminant references
- Generates detailed alignment logs and statistics
- Creates consensus sequences from cleaned alignments
- Supports different genetic codes through translation table selection
- Produces comprehensive reports of filtering decisions

## Requirements

- Python 3.6+
- Required Python packages:
  - Biopython
  - argparse
  - typing

## Installation

1. Clone this repository:
```bash
git clone https://github.com/yourusername/fasta-cleaner-contam.git
cd fasta-cleaner-contam
```

2. Install required packages:
```bash
pip install biopython
```

## Usage

### Basic Command

```bash
python fasta_cleaner_contam.py -i input_directory -o output_directory -t targets_directory -n contaminants.fasta
```

### Required Arguments

- `-i, --input_dir`: Directory containing input FASTA alignment files
- `-o, --output_dir`: Directory where output files will be saved
- `-t, --targets_dir`: Directory containing target sequence FASTA files
- `-n, --nontarget`: Path to FASTA file containing non-target contaminant sequences

### Optional Arguments

- `-c, --consensus_threshold`: Consensus threshold (0.0-1.0) for generating consensus sequences (default: 0.7)
- `--translation_table`: Genetic code table number for translation (default: 1 Standard)
  - Common options:
    - 1: Standard
    - 2: Vertebrate Mitochondrial
    - 5: Invertebrate Mitochondrial
    - 11: Bacterial and Plant Plastid

## Input File Requirements

### Input Alignments
- Must be in FASTA format
- File extensions: .fasta, .fas, or .fa
- Should contain aligned DNA sequences

### Target Sequences
- One FASTA file per input alignment
- Filename should match the input alignment (without '_align' suffix)
- Should contain a single protein sequence

### Contaminant Sequences
- Single FASTA file containing potential contaminant protein sequences
- Used for comparison against all input sequences

## Output Files

For each input FASTA file, the tool generates:

1. `*_cleaned.fasta`: Sequences that passed contamination filtering
2. `*_removed.fasta`: Sequences identified as potential contaminants
3. `*_consensus.fasta`: Consensus sequence generated from cleaned sequences
4. `*_report.csv`: Detailed report of filtering decisions including:
   - Sequence names
   - Keep/remove status
   - Target alignment scores
   - Contaminant alignment scores
   - Best matching contaminant
5. `*_detailed.log`: Comprehensive processing log including:
   - Frame translations
   - Alignment details
   - Scoring information
   - Consensus generation process

## Processing Overview

1. DNA sequences are extracted from input alignments (gaps removed)
2. Each sequence is translated in all three reading frames
3. Protein translations are aligned against:
   - Target sequence (expected protein)
   - Contaminant sequences (potential contaminants)
4. Best-scoring frame is selected for each comparison
5. Sequences are kept if they align better to target than any contaminant
6. Consensus sequence is generated from kept sequences
7. Detailed reports and logs are generated

## Error Handling

The tool includes robust error handling for:
- Missing or invalid input files
- Empty FASTA files
- Invalid translation tables
- Consensus generation issues
- File writing errors

Each error is logged appropriately, and the tool continues processing other files when possible.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.
