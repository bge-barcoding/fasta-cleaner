# FASTA Sequence Cleaner

A comprehensive Python tool for cleaning and analyzing FASTA sequence alignments using multiple filtering approaches. This tool is designed to help researchers identify and remove problematic sequences from their alignments while maintaining data quality and integrity.

## Features

The FASTA Sequence Cleaner implements a sophisticated filtering pipeline that includes:

1. **Human COX1 Contamination Detection**
   - Identifies sequences with high similarity to human COX1
   - Configurable similarity threshold
   - Helps prevent contamination from human DNA

2. **AT Content Analysis**
   - Compares AT content between sequences and consensus
   - Identifies sequences with divergent nucleotide composition
   - Customizable difference threshold

3. **Statistical Outlier Detection**
   - Uses both weighted and unweighted deviation scores
   - Position-specific residue frequency analysis
   - Adjustable percentile threshold for outlier detection

4. **Reference Sequence Comparison**
   - Optional comparison against known reference sequences
   - Supports multiple reference sequence files
   - Additional metrics for reference-based filtering

## Installation

### Prerequisites

- Python 3.6 or higher
- BioPython
- NumPy

```bash
pip install biopython numpy
```

### Installation Steps

1. Clone this repository:
```bash
git clone https://github.com/bge-barcoding/fasta-cleaner.git
cd fasta-cleaner
```

2. Install dependencies:
```bash
pip install biopython numpy
```

## Usage

### Basic Usage

```bash
python fasta_cleaner_combined.py -i input_dir -o output_dir
```

### Advanced Usage

```bash
python fasta_cleaner_combined.py \
    -i input_dir \
    -o output_dir \
    -r reference_dir \
    --human_threshold 0.95 \
    --at_difference 0.1 \
    --percentile_threshold 90.0 \
    --consensus_threshold 0.5
```

### Command Line Arguments

| Argument | Description | Default |
|----------|-------------|---------|
| `-i`, `--input_dir` | Directory containing input FASTA files | Required |
| `-o`, `--output_dir` | Output directory for processed files | Required |
| `-r`, `--reference_dir` | Directory containing reference sequences | Optional |
| `-u`, `--human_threshold` | Human COX1 similarity threshold (0-1) | 0.95 |
| `-d`, `--at_difference` | Maximum allowed AT content difference | 0.1 |
| `-p`, `--percentile_threshold` | Percentile for outlier detection (0-100) | 90.0 |
| `-c`, `--consensus_threshold` | Consensus sequence generation threshold | 0.5 |

### Filter Control Flags

| Flag | Description |
|------|-------------|
| `--disable_human` | Disable human COX1 similarity filtering |
| `--disable_at` | Disable AT content difference filtering |
| `--disable_outliers` | Disable statistical outlier detection |

## Output Files

For each input FASTA file, the tool generates:

1. `*_cleaned.fasta`: Sequences that passed all filters
2. `*_removed_all.fasta`: All removed sequences combined
3. `*_removed_human.fasta`: Sequences removed due to human similarity
4. `*_removed_at.fasta`: Sequences removed due to AT content
5. `*_removed_outlier.fasta`: Statistical outliers
6. `*_removed_reference.fasta`: Reference-based outliers
7. `*_consensus.fasta`: Generated consensus sequence
8. `*_metrics.csv`: Detailed metrics for all sequences
9. `*_log.txt`: Processing log with parameters and statistics

## Metrics and Analysis

The tool calculates comprehensive metrics for each sequence:

- Sequence length and composition
- AT content and deviation from consensus
- Human COX1 similarity scores
- Position-specific conservation scores
- Statistical deviation measures
- Reference-based metrics (if enabled)

All metrics are saved in the CSV report for further analysis.

## Example

```python
# Process a directory of FASTA files with custom thresholds
python fasta_cleaner_combined.py \
    -i /path/to/fasta/files \
    -o /path/to/output \
    -r /path/to/references \
    --human_threshold 0.90 \
    --at_difference 0.15
```

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Citation

If you use this tool in your research, please cite:

```bibtex
@software{fasta_sequence_cleaner,
  author = {Ben Price AND Daniel Parsons AND Jordan Beasley},
  title = {FASTA Sequence Cleaner},
  year = {2024},
  url = {https://github.com/bge-barcoding/fasta-cleaner}
}
```

## Acknowledgments

- Uses BioPython for sequence analysis
- Implements methods inspired by various sequence quality control approaches
- Developed to address common contamination and quality issues in sequence data

## Support

For bugs, feature requests, or questions, please open an issue on GitHub.
