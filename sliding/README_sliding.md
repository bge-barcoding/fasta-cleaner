# FASTA Cleaner with human detection and sliding AT ratios 

This Python script processes FASTA sequence alignments (e.g. from https://github.com/cmayer/MitoGeneExtractor) to identify and remove potentially contaminated or divergent sequences. 
It uses two main filtering criteria: 
- similarity to human COX1 sequences to remove human contamination
- deviation from the consensus AT content (in the local region) to remove sequences with unusual base composition (e.g. Fungi).

This builds on the script with the fixed AT ratio by allowing the consensus AT ratio to change over the length of the sequence.

## Features

The script provides several key features for sequence analysis and filtering:

- Processes multiple FASTA alignment files in batch
- Generates consensus sequences from alignments
- Filters sequences based on AT content deviation from consensus
- Removes potential human contamination by comparing to COX1
- Produces detailed reports of sequence metrics
- Preserves both kept and removed sequences for verification
- Outputs single-line FASTA format for compatibility

## Prerequisites

The script requires Python 3.7 or later and the following Python packages:

- Biopython (for sequence analysis and alignment handling)
- argparse (included in Python standard library)
- csv (included in Python standard library)
- os (included in Python standard library)

Install required packages using pip:

```bash
pip install biopython
```

## Usage

Basic usage:

```bash
python fasta_cleaner.py -i input_directory -o output_directory
```

Full usage with all options:

```bash
python fasta_cleaner.py -i input_directory -o output_directory -d 0.1 -c 0.7 -u 0.95
```

### Command Line Arguments

- `-i, --input_dir`: Directory containing input FASTA files (required)
- `-o, --output_dir`: Directory where output files will be saved (required)
- `-d, --at_difference`: Maximum allowed AT content difference from consensus (default: 0.1)
- `-c, --consensus_threshold`: Minimum frequency for consensus base calling (default: 0.7)
- `-u, --human_threshold`: Human COX1 similarity threshold for removal (default: 0.95)

### Filtering Process

The script processes each FASTA file through the following steps:

1. Reads the input alignment
2. Generates an initial consensus sequence
3. For each sequence in the alignment:
   - Calculates AT content where it overlaps with the consensus
   - Compares sequence to human COX1
   - Removes sequences that either:
     - Have AT content differing from consensus by more than the threshold
     - Show high similarity to human COX1
4. Generates a new consensus from kept sequences
5. Outputs results and detailed metrics

### Output Files

For each input FASTA file, the script generates:

- `*_cleaned.fasta`: Sequences that passed all filters
- `*_removed.fasta`: Sequences that were filtered out
- `*_cleaned_consensus.fasta`: Consensus sequence from kept sequences
- `*_report.csv`: Detailed metrics for all sequences

The CSV report includes:
- Sequence name
- Sequence length (excluding gaps)
- Sequence AT ratio
- Consensus AT ratio
- AT ratio difference
- Human similarity score
- Sequence status (kept/removed)

## Example

```bash
python fasta_cleaner.py -i sequences/raw -o sequences/cleaned -d 0.15 -c 0.6 -u 0.9
```

This command will:
- Process all FASTA files in 'sequences/raw'
- Allow sequences to differ from consensus AT content by up to 15%
- Use 60% frequency threshold for consensus calling
- Remove sequences with ≥90% similarity to human COX1
- Save results in 'sequences/cleaned'

## Interpreting Results

The script provides progress updates as it processes files:

```
Processed example.fasta:
  - Original sequences: 100
  - Kept sequences: 85
  - Removed sequences: 15
  - Saved cleaned alignment to: example_cleaned.fasta
  - Saved removed sequences to: example_removed.fasta
  - Saved consensus to: example_cleaned_consensus.fasta
  - Saved sequence report to: example_report.csv
```

Check the CSV report for detailed metrics on why specific sequences were removed.

## Notes and Recommendations

- Always examine removed sequences to ensure filtering is appropriate for your data
- Adjust thresholds based on your expected sequence diversity
- The AT difference threshold should reflect natural variation in your sequences
- The consensus threshold may need adjustment for highly variable alignments
- Keep the human threshold high (>0.9) to avoid removing non-human sequences

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.
