# FASTA Cleaner

A Python tool for cleaning and analyzing FASTA sequence alignments (e.g. from https://github.com/cmayer/MitoGeneExtractor) using multiple filtering approaches. This tool is designed to help researchers identify and remove problematic sequences from their alignments while maintaining data quality and integrity.

## Features

The FASTA Sequence Cleaner implements a sophisticated filtering pipeline that includes:

1. **Human COX1 Contamination Detection**
   - Identifies sequences with high similarity to human COX1
   - Configurable similarity threshold
   - Uses efficient local alignment for comparison
   - Helps prevent contamination from human DNA

2. **AT Content Analysis**
   - Compares AT content between sequences and consensus
   - Identifies sequences with divergent nucleotide composition
   - Supports multiple filtering modes (absolute, higher, lower)
   - Customizable difference threshold
   - Only considers overlapping regions for comparison

3. **Statistical Outlier Detection**
   - Uses both weighted and unweighted deviation scores
   - Position-specific residue frequency analysis
   - Conservation-weighted sequence comparison
   - Adjustable percentile threshold for outlier detection
   - Robust handling of gaps and missing data

4. **Reference Sequence Comparison**
   - Optional comparison against known reference sequences
   - Supports multiple reference sequence files
   - Additional metrics for reference-based filtering
   - Weighted deviation scoring based on conservation

## Installation

### Prerequisites

- Python 3.6 or higher
- BioPython
- NumPy
- typing (for type hints)

```bash
pip install biopython numpy typing
```

### Installation Steps

1. Clone this repository:
```bash
git clone https://github.com/bge-barcoding/fasta-cleaner.git
cd fasta-cleaner
```

2. Install dependencies:
```bash
pip install biopython numpy typing
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
    --at_mode absolute \
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
| `-m`, `--at_mode` | AT content filtering mode (absolute/higher/lower) | absolute |
| `-p`, `--percentile_threshold` | Percentile for outlier detection (0-100) | 90.0 |
| `-c`, `--consensus_threshold` | Consensus sequence generation threshold | 0.5 |

### AT Content Filtering Modes

The tool supports three modes for AT content filtering:

- `absolute`: Removes sequences if AT content differs from consensus by more than threshold in either direction
- `higher`: Removes only sequences with AT content above consensus + threshold (i.e. AT is too high)
- `lower`: Removes only sequences with AT content below consensus - threshold (i.e. AT is too low)

### Filter Control Flags

| Flag | Description |
|------|-------------|
| `--disable_human` | Disable human COX1 similarity filtering |
| `--disable_at` | Disable AT content difference filtering |
| `--disable_outliers` | Disable statistical outlier detection |

## Output Files

For each input FASTA file, the tool generates:

1. `*_cleaned.fasta`: Sequences that passed all filters, ordered by start position
2. `*_removed_all.fasta`: All removed sequences combined into one file
3. `*_removed_human.fasta`: Sequences removed due to human similarity
4. `*_removed_at.fasta`: Sequences removed due to AT content
5. `*_removed_outlier.fasta`: Sequences removed as statistical outliers
6. `*_removed_reference.fasta`:  Sequences removed as reference-based outliers
7. `*_consensus.fasta`: Final consensus sequence
8. `*_metrics.csv`: Detailed metrics for all sequences
9. `*_log.txt`: Processing log with parameters and statistics
10. `*_ordered_annotated.fasta`: All original sequences with fate annotations, ordered by start position

## Metrics and Analysis

The tool calculates comprehensive metrics for each sequence:

- Sequence length and composition
- AT content and deviation from consensus
- Human COX1 similarity scores using local alignment
- Position-specific conservation scores
- Weighted and unweighted deviation measures
- Conservation-based statistical scores
- Reference-based metrics (if enabled)
- Gap handling and position-specific frequencies

All metrics are saved in the CSV report for further analysis.

## Sequence Processing Pipeline

The filtering pipeline processes sequences in this specific order:

1. Remove sequences with high human COX1 similarity
2. Filter sequences with divergent AT content
3. Remove statistical outliers
4. Compare against reference sequences (if provided)

After each filtering step:
- A new consensus sequence is generated from remaining sequences
- New position-specific frequencies are calculated
- New metrics are computed for all remaining sequences

## Example

```python
# Process a directory of FASTA files with custom thresholds and AT mode lower
python fasta_cleaner_combined.py \
    -i /path/to/fasta/files \
    -o /path/to/output \
    -r /path/to/references \
    --human_threshold 0.90 \
    --at_difference 0.15 \
    --at_mode lower \
    --percentile_threshold 95.0
```

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request. 

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Citation

If you use this tool in your research, please cite:

```bibtex
@software{fasta_cleaner,
  author = {Ben Price AND Daniel Parsons AND Jordan Beasley AND Claude Sonnet},
  title = {FASTA Sequence Cleaner},
  version = {1.0.0},
  year = {2024},
  url = {https://github.com/bge-barcoding/fasta-cleaner},
  note = {Implements multiple sequence filtering approaches with position-specific analysis}
}
```

## Acknowledgments

- Uses BioPython for sequence analysis
- Implements methods inspired by various sequence quality control approaches
- Developed to address common contamination and quality issues in sequence data

## Support

For bugs, feature requests, or questions, please open an issue on GitHub.
