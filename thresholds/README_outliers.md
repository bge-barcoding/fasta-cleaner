# FASTA Cleaner with outlier detection

A robust Python tool for identifying and filtering outlier sequences in FASTA alignment files using consensus-based thresholding. This tool is particularly useful for cleaning multiple sequence alignments by removing sequences that significantly deviate from the consensus.

## Features

The script provides comprehensive sequence analysis capabilities including:

- Consensus sequence generation with customizable threshold
- Position-specific residue frequency calculation
- Multiple scoring methods for sequence deviation:
  - Unweighted deviation from consensus
  - Weighted deviation considering position-specific conservation
  - Optional comparison against reference sequences
- Automatic outlier detection using statistical thresholding
- Detailed logging of the analysis process
- Generation of cleaned and filtered sequence sets

## Requirements

- Python 3.6 or higher
- Required Python packages:
  - BioPython
  - NumPy
  - argparse
  - typing

## Installation

1. Clone this repository:
```bash
git clone https://github.com/yourusername/fasta-cleaner.git
cd fasta-cleaner
```

2. Install required packages:
```bash
pip install biopython numpy
```

## Usage

The script can be run from the command line with several options:

```bash
python identify_outliers_thresholding.py -i input_directory -o output_directory [-r reference_directory] [-c consensus_threshold]
```

### Required Arguments

- `-i, --input_dir`: Directory containing input FASTA alignment files
- `-o, --output_dir`: Directory where output files will be saved

### Optional Arguments

- `-r, --reference_dir`: Directory containing reference FASTA files (must have "_reference.fasta" suffix)
- `-c, --consensus_threshold`: Threshold for consensus sequence generation (default: 0.7)

## Output Files

For each input FASTA file, the script generates several output files:

1. `*_scores.csv`: Contains deviation scores and fate (kept/removed) for each sequence
2. `*_cleaned.fasta`: FASTA file containing sequences that passed the filtering
3. `*_removed.fasta`: FASTA file containing sequences that were filtered out
4. `*_consensus.fasta`: Generated consensus sequence
5. `*_log.txt`: Detailed log file containing processing information and statistics

## How It Works

1. **Sequence Analysis**:
   - Calculates position-specific residue frequencies
   - Generates consensus sequence using specified threshold
   - Computes deviation scores for each sequence

2. **Outlier Detection**:
   - Calculates statistical thresholds (75th percentile) from unweighted deviation scores
   - Classifies sequences as outliers if they exceed the threshold

3. **Scoring Methods**:
   - Unweighted Deviation: Simple proportion of mismatches
   - Weighted Deviation: Considers position-specific conservation weights
   - Reference-based Scoring: Optional comparison against reference sequences

## Example

```bash
python identify_outliers_thresholding.py \
  -i /path/to/alignments \
  -o /path/to/output \
  -r /path/to/references \
  -c 0.8
```

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.
