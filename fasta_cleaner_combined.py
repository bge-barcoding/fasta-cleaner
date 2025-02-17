import os
import sys
import csv
import argparse
from datetime import datetime
from typing import List, Dict, Tuple, Optional, Any, NamedTuple
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO
from Bio.Align import AlignInfo
import numpy as np
from collections import Counter
from collections import defaultdict
import traceback

"""FASTA Sequence Cleaner and Analysis Tool

A comprehensive tool for cleaning and analysing FASTA sequence alignments using multiple filtering approaches:
- Human COX1 similarity detection
- AT-content analysis
- Statistical outlier detection
- Reference sequence comparison

The tool processes each alignment file through a sequential filtering pipeline:
1. Remove sequences with high human COX1 similarity
2. Filter sequences with divergent AT content
3. Remove statistical outliers
4. Compare against reference sequences (if provided)

Input:
------
Required:
    - Directory containing FASTA alignment files (.fasta, .fas, or .fa)
    - Output directory for results

Optional:
    - Directory containing reference sequences (named as original_filename_reference.fasta)

Output:
-------
For each input file, generates:
    - {basename}_cleaned.fasta: Filtered sequences that passed all criteria
    - {basename}_consensus.fasta: Consensus sequence generated from cleaned alignment
    - {basename}_metrics.csv: Detailed metrics for all sequences
    - {basename}_ordered_annotated.fasta: All sequences with filtering annotations
    - {basename}_removed_all.fasta: All removed sequences combined
    - {basename}_removed_human.fasta: Sequences removed due to human similarity
    - {basename}_removed_at.fasta: Sequences removed due to AT content
    - {basename}_removed_outlier.fasta: Sequences removed as statistical outliers
    - {basename}_removed_reference.fasta: Sequences removed based on reference comparison
    - {basename}_log.txt: Detailed processing log

Options:
--------
Filter Enablement:
    --disable_human: Skip human COX1 similarity filtering
    --disable_at: Skip AT content filtering
    --disable_outliers: Skip statistical outlier detection

Human COX1 Filtering:
    -u, --human_threshold: Similarity threshold (0.0-1.0, default: 0.95)
        Sequences with similarity >= threshold are removed

AT Content Filtering:
    -d, --at_difference: Maximum allowed AT content difference (0.0-1.0, default: 0.1)
    -m, --at_mode: Filtering mode (default: 'absolute')
        - 'absolute': Remove if AT content differs from consensus by more than threshold
        - 'higher': Remove only sequences with AT content above consensus + threshold
        - 'lower': Remove only sequences with AT content below consensus - threshold

Statistical Outlier Detection:
    -p, --percentile_threshold: Percentile for outlier detection (0.0-100.0, default: 90.0)
        Sequences with deviation scores above this percentile are removed

Consensus Generation:
    -c, --consensus_threshold: Threshold for consensus sequence generation (0.0-1.0, default: 0.5)
        Minimum frequency required to call a consensus base

Reference Sequence Comparison:
    -r, --reference_dir: Directory containing reference sequences
        Reference files should be named as {original_filename}_reference.fasta

Usage:
------
Basic usage:
    python fasta_cleaner_combined.py -i input_dir -o output_dir

With all options:
    python fasta_cleaner_combined.py -i input_dir -o output_dir \\
        -r reference_dir \\
        -u 0.95 \\
        -d 0.1 \\
        -m absolute \\
        -p 90.0 \\
        -c 0.5

Disable specific filters:
    python fasta_cleaner_combined.py -i input_dir -o output_dir \\
        --disable_human --disable_at

Requirements:
------------
Python packages:
    - BioPython
    - NumPy
    - typing
    - argparse

Example:
--------
python fasta_cleaner_combined.py \\
    -i /path/to/fasta/files \\
    -o /path/to/output \\
    -r /path/to/references \\
    -u 0.98 \\
    -d 0.15 \\
    -m absolute \\
    -p 95.0
"""

# Reference human COX1 sequence
HUMAN_COX1 = """ATGTTCGCCGACCGTTGACTATTCTCTACAAACCACAAAGACATTGGAACACTATACCTATTATTCGGCG
CATGAGCTGGAGTCCTAGGCACAGCTCTAAGCCTCCTTATTCGAGCCGAGCTGGGCCAGCCAGGCAACCT
TCTAGGTAACGACCACATCTACAACGTTATCGTCACAGCCCATGCATTTGTAATAATCTTCTTCATAGTA
ATACCCATCATAATCGGAGGCTTTGGCAACTGACTAGTTCCCCTAATAATCGGTGCCCCCGATATGGCGT
TTCCCCGCATAAACAACATAAGCTTCTGACTCTTACCTCCCTCTCTCCTACTCCTGCTCGCATCTGCTAT
AGTGGAGGCCGGAGCAGGAACAGGTTGAACAGTCTACCCTCCCTTAGCAGGGAACTACTCCCACCCTGGA
GCCTCCGTAGACCTAACCATCTTCTCCTTACACCTAGCAGGTGTCTCCTCTATCTTAGGGGCCATCAATT
TCATCACAACAATTATCAATATAAAACCCCCTGCCATAACCCAATACCAAACGCCCCTCTTCGTCTGATC
CGTCCTAATCACAGCAGTCCTACTTCTCCTATCTCTCCCAGTCCTAGCTGCTGGCATCACTATACTACTA
ACAGACCGCAACCTCAACACCACCTTCTTCGACCCCGCCGGAGGAGGAGACCCCATTCTATACCAACACC
TATTCTGATTTTTCGGTCACCCTGAAGTTTATATTCTTATCCTACCAGGCTTCGGAATAATCTCCCATAT
TGTAACTTACTACTCCGGAAAAAAAGAACCATTTGGATACATAGGTATGGTCTGAGCTATGATATCAATT
GGCTTCCTAGGGTTTATCGTGTGAGCACACCATATATTTACAGTAGGAATAGACGTAGACACACGAGCAT
ATTTCACCTCCGCTACCATAATCATCGCTATCCCCACCGGCGTCAAAGTATTTAGCTGACTCGCCACACT
CCACGGAAGCAATATGAAATGATCTGCTGCAGTGCTCTGAGCCCTAGGATTCATCTTTCTTTTCACCGTA
GGTGGCCTGACTGGCATTGTATTAGCAAACTCATCACTAGACATCGTACTACACGACACGTACTACGTTG
TAGCCCACTTCCACTATGTCCTATCAATAGGAGCTGTATTTGCCATCATAGGAGGCTTCATTCACTGATT
TCCCCTATTCTCAGGCTACACCCTAGACCAAACCTACGCCAAAATCCATTTCACTATCATATTCATCGGC
GTAAATCTAACTTTCTTCCCACAACACTTTCTCGGCCTATCCGGAATGCCCCGACGTTACTCGGACTACC
CCGATGCATACACCACATGAAACATCCTATCATCTGTAGGCTCATTCATTTCTCTAACAGCAGTAATATT
AATAATTTTCATGATTTGAGAAGCCTTCGCTTCGAAGCGAAAAGTCCTAATAGTAGAAGAACCCTCCATA
AACCTGGAGTGACTATATGGATGCCCCCCACCCTACCACACATTCGAAGAACCCGTATACATAAAATCTA
GA""".replace("\n", "")

def parse_arguments() -> argparse.Namespace:
    """
    Set up command line argument parsing with comprehensive options for all filtering methods.
    Includes validation of input parameters and directory paths.
    """
    parser = argparse.ArgumentParser(
        description='Advanced FASTA sequence cleaner with multiple filtering methods: '
                   'human COX1 similarity, AT content, and statistical outlier detection.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Required arguments
    parser.add_argument(
        '-i', '--input_dir',
        required=True,
        help='Directory containing input FASTA alignment files'
    )
    
    parser.add_argument(
        '-o', '--output_dir',
        required=True,
        help='Directory where output files will be saved'
    )
    
    # Filter enablement flags
    parser.add_argument(
        '--disable_human',
        action='store_true',
        help='Disable human COX1 similarity filtering'
    )
    
    parser.add_argument(
        '--disable_at',
        action='store_true',
        help='Disable AT content difference filtering'
    )
    
    parser.add_argument(
        '--disable_outliers',
        action='store_true',
        help='Disable statistical outlier detection'
    )
    
    # Human COX1 filtering parameters
    parser.add_argument(
        '-u', '--human_threshold',
        type=float,
        default=0.95,
        help='Human COX1 similarity threshold (0.0-1.0). Sequences with similarity >= threshold are removed'
    )
    
    # AT content filtering parameters
    parser.add_argument(
        '-d', '--at_difference',
        type=float,
        default=0.1,
        help='Maximum allowed AT content difference from consensus (0.0-1.0)'
    )
    
    parser.add_argument(
        '-m','--at_mode',
        type=str,
        choices=['absolute', 'higher', 'lower'],
        default='absolute',
        help='AT content filtering mode: "absolute" removes sequences if AT content differs '
             'from consensus by more than threshold in either direction, "higher" removes only '
             'sequences with AT content above consensus + threshold, "lower" removes only '
             'sequences with AT content below consensus - threshold'
    )
    
    # Statistical outlier detection parameters
    parser.add_argument(
        '-p', '--percentile_threshold',
        type=float,
        default=90.0,
        help='Percentile threshold for statistical outlier detection (0.0-100.0)'
    )
    
    # Consensus generation parameters
    parser.add_argument(
        '-c', '--consensus_threshold',
        type=float,
        default=0.5,
        help='Threshold for consensus sequence generation (0.0-1.0)'
    )
    
    # Optional reference sequence directory
    parser.add_argument(
        '-r', '--reference_dir',
        help='Optional directory containing reference FASTA files (named same as input files but with "_reference.fasta" suffix)'
    )
    
    args = parser.parse_args()
    
    # Validate directory paths
    if not os.path.isdir(args.input_dir):
        parser.error(f"Input directory does not exist: {args.input_dir}")
    
    if args.reference_dir and not os.path.isdir(args.reference_dir):
        parser.error(f"Reference directory does not exist: {args.reference_dir}")
    
    # Validate thresholds
    if not 0 <= args.human_threshold <= 1:
        parser.error("Human threshold must be between 0.0 and 1.0")
    
    if not 0 <= args.at_difference <= 1:
        parser.error("AT difference threshold must be between 0.0 and 1.0")
    
    if not 0 <= args.percentile_threshold <= 100:
        parser.error("Percentile threshold must be between 0.0 and 100.0")
    
    if not 0 < args.consensus_threshold <= 1:
        parser.error("Consensus threshold must be between 0.0 and 1.0")
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    return args

# Custom type definitions for improved code clarity
SequenceRecord = SeqRecord
SequenceAlignment = MultipleSeqAlignment

def get_sequence_start_position(sequence: str) -> int:
    """
    Find the position of the first non-gap character in a sequence.
    
    Args:
        sequence: Aligned sequence string potentially containing gaps
        
    Returns:
        int: Position of first non-gap character (0-based) or length of sequence if all gaps
    """
    for i, char in enumerate(sequence):
        if char != '-':
            return i
    return len(sequence)

def sort_sequences_by_position(records: List[SeqRecord]) -> List[SeqRecord]:
    """
    Sort sequence records by their starting position in the alignment.
    
    Args:
        records: List of SeqRecord objects to sort
        
    Returns:
        List[SeqRecord]: Sorted list of records
    """
    # Create tuples of (start_position, original_index, record) to maintain stable sort
    indexed_records = [
        (get_sequence_start_position(str(record.seq)), i, record)
        for i, record in enumerate(records)
    ]
    
    # Sort by start position, using original index as tiebreaker
    sorted_records = [
        record for _, _, record in sorted(indexed_records)
    ]
    
    return sorted_records

def create_annotated_sequence_record(record: SeqRecord, fate: str) -> SeqRecord:
    """
    Create a new sequence record with the fate prefix added to its ID.
    
    Args:
        record: Original SeqRecord object
        fate: String indicating the sequence's fate (e.g., 'kept', 'removed_human')
        
    Returns:
        SeqRecord: New record with modified ID
    """
    # Create a deep copy to avoid modifying the original
    new_record = SeqRecord(
        seq=record.seq,
        id=f"{fate}_{record.id}",
        name=record.name,
        description=record.description
    )
    return new_record

def write_ordered_annotated_alignment(
    kept_records: List[SeqRecord],
    removed_records: Dict[str, List[SeqRecord]],
    output_path: str
) -> None:
    """
    Write all sequences to a FASTA file with fate annotations and ordered by start position.
    
    Args:
        kept_records: List of sequences that passed all filters
        removed_records: Dictionary mapping removal reasons to lists of removed sequences
        output_path: Path to write the output file
    """
    # Create annotated records for all sequences
    all_records = []
    
    # Add kept sequences with 'kept' prefix
    all_records.extend(
        create_annotated_sequence_record(record, 'kept')
        for record in kept_records
    )
    
    # Add removed sequences with appropriate prefixes
    for reason, records in removed_records.items():
        prefix = f"removed_{reason}"
        all_records.extend(
            create_annotated_sequence_record(record, prefix)
            for record in records
        )
    
    # Sort all records by start position
    sorted_records = sort_sequences_by_position(all_records)
    
    # Write to file
    write_sequences_to_fasta(sorted_records, output_path)

def find_best_alignment(query: str, reference: str, min_overlap: int = 20) -> tuple[int, int, int]:
    """
    Find the best local alignment between a query sequence and a reference sequence.
    
    Args:
        query: Query sequence string (without gaps)
        reference: Reference sequence string
        min_overlap: Minimum required overlap length (default: 20)
    
    Returns:
        Tuple of (start_pos, end_pos, match_count):
            start_pos: Best starting position in reference
            end_pos: Best ending position in reference
            match_count: Number of matches in best alignment
    """
    query_len = len(query)
    ref_len = len(reference)
    
    best_start = 0
    best_end = 0
    max_matches = 0
    
    # Try all possible alignments of query against reference
    for start in range(ref_len - min_overlap + 1):
        # Don't look past end of reference
        end = min(start + query_len, ref_len)
        
        # Count matches in this alignment window
        ref_segment = reference[start:end]
        query_segment = query[:end-start]
        
        matches = sum(1 for q, r in zip(query_segment, ref_segment) 
                     if q == r and q != '-' and r != '-')
        
        # Update best alignment if we found more matches
        if matches > max_matches:
            max_matches = matches
            best_start = start
            best_end = end

    return best_start, best_end, max_matches

def calculate_sequence_similarity(seq1: str, seq2: str, min_overlap: int = 20) -> float:
    """
    Calculate sequence similarity between two sequences by finding the best local alignment.
    The query sequence (seq1) is stripped of gaps and aligned against the reference (seq2).
    Similarity is calculated only over the best matching aligned region.
    
    Args:
        seq1: Query sequence string (can contain gaps)
        seq2: Reference sequence string (can contain gaps)
        min_overlap: Minimum required overlap length (default: 20)
    
    Returns:
        float: Similarity ratio between 0 and 1 for best matching region
               Returns 0.0 if no valid alignment found
    """
    # Remove gaps and convert to uppercase
    query = seq1.replace('-', '').upper()
    reference = seq2.replace('-', '').upper()
    
    # Check sequence lengths
    if len(query) < min_overlap or len(reference) < min_overlap:
        return 0.0
        
    # Find best local alignment
    start, end, matches = find_best_alignment(query, reference, min_overlap)
    
    # Calculate similarity over aligned region
    aligned_length = end - start
    if aligned_length < min_overlap:
        return 0.0
        
    return matches / aligned_length


    """
    Calculate sequence similarity between two sequences, ignoring gaps.
    Uses a position-by-position comparison of ungapped regions.
    
    Args:
        seq1: First sequence string
        seq2: Second sequence string
    
    Returns:
        float: Similarity ratio between 0 and 1
    """
    # Remove gaps and convert to uppercase for comparison
    seq1_nogaps = seq1.replace('-', '').upper()
    seq2_nogaps = seq2.replace('-', '').upper()
    
    # Get the length of the shorter sequence
    min_length = min(len(seq1_nogaps), len(seq2_nogaps))
    
    if min_length == 0:
        return 0.0
    
    # Count matches in the overlapping region
    matches = sum(1 for i in range(min_length) 
                 if seq1_nogaps[i] == seq2_nogaps[i])
    
    return matches / min_length

def calculate_at_content(sequence: str) -> float:
    """
    Calculate AT content for a sequence, ignoring gaps.
    
    Args:
        sequence: Input sequence string
    
    Returns:
        float: AT content ratio between 0 and 1
    """
    sequence_no_gaps = sequence.replace('-', '').upper()
    
    if not sequence_no_gaps:
        return 0.0
    
    at_count = sequence_no_gaps.count('A') + sequence_no_gaps.count('T')
    return at_count / len(sequence_no_gaps)

def compare_at_content(consensus_seq: str, query_seq: str) -> Tuple[float, float, float]:
    """
    Compare AT content between consensus and query sequence in overlapping regions.
    Only considers positions where both sequences have bases (not gaps).
    
    Args:
        consensus_seq: Consensus sequence string
        query_seq: Query sequence string
    
    Returns:
        Tuple containing:
        - query AT content
        - consensus AT content
        - absolute difference between AT contents
    """
    consensus_seq = consensus_seq.upper()
    query_seq = query_seq.upper()
    
    # Find positions where both sequences have bases
    overlapping_positions = [
        (c, q) for c, q in zip(consensus_seq, query_seq)
        if c != '-' and q != '-'
    ]
    
    if not overlapping_positions:
        return 0.0, 0.0, 0.0
    
    # Split into consensus and query sequences
    consensus_bases = ''.join(c for c, _ in overlapping_positions)
    query_bases = ''.join(q for _, q in overlapping_positions)
    
    # Calculate AT content for overlapping regions
    consensus_at = calculate_at_content(consensus_bases)
    query_at = calculate_at_content(query_bases)
    
    return query_at, consensus_at, abs(query_at - consensus_at)

def calculate_position_frequencies(sequences: List[str]) -> List[Dict[str, float]]:
    """
    Calculate residue frequencies at each position in the alignment.
    Gaps are excluded from frequency calculations.
    
    Args:
        sequences: List of aligned sequence strings
    
    Returns:
        List of dictionaries containing residue frequencies for each position
    """
    if not sequences:
        raise ValueError("No sequences provided")
    
    align_length = len(sequences[0])
    position_freqs = []
    
    for i in range(align_length):
        # Get residues at this position, excluding gaps
        residues = [seq[i] for seq in sequences if seq[i] != '-']
        
        if residues:
            counts = Counter(residues)
            total = sum(counts.values())
            frequencies = {res: count/total for res, count in counts.items()}
        else:
            frequencies = {}
            
        position_freqs.append(frequencies)
    
    return position_freqs

def calculate_weighted_deviation(sequence: str, reference: str, 
                              frequencies: List[Dict[str, float]]) -> float:
    """
    Calculate weighted deviation score between sequence and reference.
    Weights are based on position-specific conservation.
    
    Args:
        sequence: Query sequence string
        reference: Reference sequence string
        frequencies: List of position-specific residue frequencies
    
    Returns:
        float: Weighted deviation score between 0 and 1
    """
    if len(sequence) != len(reference) != len(frequencies):
        raise ValueError("Sequence, reference, and frequencies must all have same length")
    
    total_score = 0.0
    total_weight = 0.0
    
    for seq_res, ref_res, pos_freqs in zip(sequence, reference, frequencies):
        if seq_res != '-' and ref_res != '-':
            # Weight is based on reference residue conservation
            conservation_weight = pos_freqs.get(ref_res, 0)
            total_weight += conservation_weight
            
            if seq_res != ref_res:
                total_score += conservation_weight
    
    return total_score / total_weight if total_weight > 0 else 0.0

def calculate_unweighted_deviation(sequence: str, reference: str) -> float:
    """
    Calculate unweighted deviation score between sequence and reference.
    Simply counts proportion of mismatches in ungapped positions.
    
    Args:
        sequence: Query sequence string
        reference: Reference sequence string
    
    Returns:
        float: Unweighted deviation score between 0 and 1
    """
    if len(sequence) != len(reference):
        raise ValueError("Sequence and reference must be same length")
    
    differences = 0
    valid_positions = 0
    
    for seq_res, ref_res in zip(sequence, reference):
        if seq_res != '-' and ref_res != '-':
            valid_positions += 1
            if seq_res != ref_res:
                differences += 1
    
    return differences / valid_positions if valid_positions > 0 else 0.0

def analyze_sequence_metrics(record: SeqRecord, 
                           consensus_seq: Optional[str] = None,
                           frequencies: Optional[List[Dict[str, float]]] = None) -> Dict[str, float]:
    """
    Calculate comprehensive sequence metrics including length, AT content,
    and deviations from consensus if provided.
    
    Args:
        record: SeqRecord object containing sequence
        consensus_seq: Optional consensus sequence string
        frequencies: Optional list of position-specific frequencies
    
    Returns:
        Dictionary containing calculated metrics
    """
    sequence = str(record.seq).upper()
    sequence_no_gaps = sequence.replace('-', '')
    
    metrics = {
        'length': len(sequence_no_gaps),
        'at_content': calculate_at_content(sequence_no_gaps),
        'human_similarity': calculate_sequence_similarity(sequence, HUMAN_COX1)
    }
    
    if consensus_seq is not None:
        # Calculate AT content comparison
        query_at, cons_at, at_diff = compare_at_content(consensus_seq, sequence)
        metrics.update({
            'consensus_at': cons_at,
            'at_difference': at_diff,
            'unweighted_deviation': calculate_unweighted_deviation(sequence, consensus_seq)
        })
        
        # Add weighted deviation if frequencies are provided
        if frequencies is not None:
            metrics['weighted_deviation'] = calculate_weighted_deviation(
                sequence, consensus_seq, frequencies
            )
    
    return metrics

class FilterResult(NamedTuple):
    """Container for sequence filtering results"""
    kept_records: List[SeqRecord]
    removed_records: Dict[str, List[SeqRecord]]  # Maps reason -> list of removed records
    metrics: Dict[str, Dict[str, float]]  # Maps sequence ID -> metric values
    consensus_seq: str  # Final consensus sequence
    frequencies: List[Dict[str, float]]  # Position-specific residue frequencies
    reference_metrics: Optional[Dict[str, float]] = None  # Optional reference-based metrics

def generate_consensus_sequence(alignment: MultipleSeqAlignment, 
                              threshold: float = 0.5) -> Tuple[str, List[Dict[str, float]]]:
    """
    Generate consensus sequence and calculate position-specific frequencies.
    
    Args:
        alignment: MultipleSeqAlignment object
        threshold: Minimum frequency threshold for consensus base calling
    
    Returns:
        Tuple containing:
        - Consensus sequence string
        - List of position-specific residue frequencies
    """
    if not alignment:
        raise ValueError("Empty alignment provided")
    
    # Extract sequences as strings
    sequences = [str(record.seq) for record in alignment]
    
    # Calculate position-specific frequencies
    frequencies = calculate_position_frequencies(sequences)
    
    # Generate consensus sequence
    consensus = []
    for pos_freqs in frequencies:
        if pos_freqs:
            # Get most common residue and its frequency
            most_common = max(pos_freqs.items(), key=lambda x: x[1])
            if most_common[1] >= threshold:
                consensus.append(most_common[0])
            else:
                consensus.append('-')
        else:
            consensus.append('-')
    
    return ''.join(consensus), frequencies

def filter_human_sequences(records: List[SeqRecord], 
                         threshold: float = 0.95) -> Tuple[List[SeqRecord], List[SeqRecord]]:
    """
    Filter out sequences with high similarity to human COX1.
    
    Args:
        records: List of SeqRecord objects
        threshold: Similarity threshold for removal
    
    Returns:
        Tuple containing:
        - List of kept records
        - List of removed records
    """
    kept = []
    removed = []
    
    for record in records:
        similarity = calculate_sequence_similarity(str(record.seq), HUMAN_COX1)
        if similarity >= threshold:
            removed.append(record)
        else:
            kept.append(record)
    
    return kept, removed

def filter_at_content(records: List[SeqRecord],
                     consensus_seq: str,
                     threshold: float = 0.1,
                     mode: str = 'absolute') -> Tuple[List[SeqRecord], List[SeqRecord]]:
    """
    Filter sequences based on AT content difference from consensus.
    
    Args:
        records: List of SeqRecord objects
        consensus_seq: Consensus sequence string
        threshold: Maximum allowed AT content difference
        mode: Filtering mode ('absolute', 'higher', or 'lower')
            - 'absolute': Remove if AT content differs from consensus by more than threshold
            - 'higher': Remove if AT content exceeds consensus by more than threshold
            - 'lower': Remove if AT content is below consensus by more than threshold
    
    Returns:
        Tuple containing:
        - List of kept records
        - List of removed records
    """
    kept = []
    removed = []
    
    for record in records:
        query_at, cons_at, at_diff = compare_at_content(consensus_seq, str(record.seq))
        # Calculate signed difference (positive means query has higher AT content)
        signed_diff = query_at - cons_at
        
        should_remove = False
        if mode == 'absolute':
            should_remove = abs(signed_diff) > threshold
        elif mode == 'higher':
            should_remove = signed_diff > threshold
        elif mode == 'lower':
            should_remove = signed_diff < -threshold
        
        if should_remove:
            removed.append(record)
        else:
            kept.append(record)
    
    return kept, removed

def filter_statistical_outliers(records: List[SeqRecord],
                              consensus_seq: str,
                              frequencies: List[Dict[str, float]],
                              percentile: float = 90.0) -> Tuple[List[SeqRecord], List[SeqRecord]]:
    """
    Filter sequences based on statistical outlier detection.
    Uses both weighted and unweighted deviation scores.
    
    Args:
        records: List of SeqRecord objects
        consensus_seq: Consensus sequence string
        frequencies: List of position-specific frequencies
        percentile: Percentile threshold for outlier detection
    
    Returns:
        Tuple containing:
        - List of kept records
        - List of removed records
    """
    # Calculate deviation scores
    unweighted_scores = []
    weighted_scores = []
    
    for record in records:
        sequence = str(record.seq)
        unweighted = calculate_unweighted_deviation(sequence, consensus_seq)
        weighted = calculate_weighted_deviation(sequence, consensus_seq, frequencies)
        
        unweighted_scores.append(unweighted)
        weighted_scores.append(weighted)
    
    # Calculate thresholds (using non-zero scores)
    unweighted_threshold = np.percentile([s for s in unweighted_scores if s > 0], percentile)
    weighted_threshold = np.percentile([s for s in weighted_scores if s > 0], percentile)
    
    # Filter sequences
    kept = []
    removed = []
    
    for record, unw_score, w_score in zip(records, unweighted_scores, weighted_scores):
        if unw_score > unweighted_threshold or w_score > weighted_threshold:
            removed.append(record)
        else:
            kept.append(record)
    
    return kept, removed

def get_reference_sequence(reference_file: str) -> str:
    """
    Extract and validate reference sequence from a FASTA file.
    
    Args:
        reference_file: Path to reference FASTA file
    
    Returns:
        Reference sequence string
    
    Raises:
        ValueError: If no sequences found or error reading file
    """
    try:
        references = list(SeqIO.parse(reference_file, "fasta"))
        if not references:
            raise ValueError("No sequences found in reference file")
        if len(references) > 1:
            print(f"Warning: Multiple sequences found in {reference_file}, using first one")
        return str(references[0].seq)
    except Exception as e:
        raise ValueError(f"Error reading reference file: {str(e)}")

def calculate_reference_metrics(sequence: str, 
                             reference_seq: str,
                             frequencies: List[Dict[str, float]]) -> Dict[str, float]:
    """
    Calculate comprehensive reference-based metrics for a sequence.
    
    Args:
        sequence: Query sequence string
        reference_seq: Reference sequence string
        frequencies: Position-specific residue frequencies
    
    Returns:
        Dictionary of reference-based metrics
    """
    return {
        'reference_unweighted_deviation': calculate_unweighted_deviation(sequence, reference_seq),
        'reference_weighted_deviation': calculate_weighted_deviation(sequence, reference_seq, frequencies),
        'reference_at_difference': abs(calculate_at_content(sequence) - calculate_at_content(reference_seq))
    }

def filter_reference_outliers(records: List[SeqRecord],
                            reference_seq: str,
                            frequencies: List[Dict[str, float]],
                            percentile: float = 90.0) -> Tuple[List[SeqRecord], List[SeqRecord], Dict[str, float]]:
    """
    Filter sequences based on deviation from reference sequence.
    
    Args:
        records: List of SeqRecord objects
        reference_seq: Reference sequence string
        frequencies: Position-specific residue frequencies
        percentile: Percentile threshold for outlier detection
    
    Returns:
        Tuple containing:
        - List of kept records
        - List of removed records
        - Dictionary of reference-based metrics
    """
    # Calculate metrics for all sequences
    sequence_metrics = {
        record.id: calculate_reference_metrics(
            str(record.seq), reference_seq, frequencies
        ) for record in records
    }
    
    # Calculate thresholds using non-zero scores
    unweighted_scores = [
        metrics['reference_unweighted_deviation'] 
        for metrics in sequence_metrics.values()
        if metrics['reference_unweighted_deviation'] > 0
    ]
    weighted_scores = [
        metrics['reference_weighted_deviation']
        for metrics in sequence_metrics.values()
        if metrics['reference_weighted_deviation'] > 0
    ]
    
    # Calculate percentile thresholds
    unweighted_threshold = np.percentile(unweighted_scores, percentile)
    weighted_threshold = np.percentile(weighted_scores, percentile)
    
    # Filter sequences
    kept = []
    removed = []
    
    for record in records:
        metrics = sequence_metrics[record.id]
        if (metrics['reference_unweighted_deviation'] > unweighted_threshold or
            metrics['reference_weighted_deviation'] > weighted_threshold):
            removed.append(record)
        else:
            kept.append(record)
    
    # Calculate aggregate metrics
    aggregate_metrics = {
        'mean_unweighted_deviation': np.mean(unweighted_scores),
        'mean_weighted_deviation': np.mean(weighted_scores),
        'unweighted_threshold': unweighted_threshold,
        'weighted_threshold': weighted_threshold
    }
    
    return kept, removed, aggregate_metrics

def apply_filters(alignment: MultipleSeqAlignment,
                 consensus_threshold: float = 0.5,
                 human_threshold: float = 0.95,
                 at_threshold: float = 0.1,
                 at_mode: str = 'absolute',
                 outlier_percentile: float = 90.0,
                 reference_seq: Optional[str] = None,
                 enable_human: bool = True,
                 enable_at: bool = True,
                 enable_outliers: bool = True,
                 enable_reference: bool = True) -> FilterResult:
    """
    Apply filtering methods in strict order with consensus recalculation between steps.
    
    The filtering pipeline follows this exact order:
    1. Human COX1 similarity (if enabled)
    2. AT content difference (if enabled)
    3. Statistical outliers (if enabled)
    4. Reference sequence comparison (if enabled and provided)
    
    After each filtering step, if sequences were removed:
    - A new consensus sequence is generated from remaining sequences
    - New position-specific frequencies are calculated
    - New metrics are computed for all remaining sequences
    
    Args:
        alignment: Input sequence alignment
        consensus_threshold: Threshold for consensus generation (0-1)
        human_threshold: Human similarity threshold for removal (0-1)
        at_threshold: Maximum allowed AT content difference (0-1)
        at_mode: AT content filtering mode ('absolute', 'higher', or 'lower')
        outlier_percentile: Percentile threshold for outlier detection (0-100)
        reference_seq: Optional reference sequence string
        enable_human: Enable human similarity filtering
        enable_at: Enable AT content filtering
        enable_outliers: Enable statistical outlier detection
        enable_reference: Enable reference sequence filtering
    
    Returns:
        FilterResult containing:
        - kept_records: Sequences that passed all filters
        - removed_records: Dict mapping removal reasons to removed sequences
        - metrics: Dict mapping sequence IDs to their metric values
        - consensus_seq: Final consensus sequence after all filtering
    """
    # Initialize our tracking structures
    current_records = list(alignment)
    removed_records = {
        'human_similar': [],
        'at_difference': [],
        'statistical_outlier': [],
        'reference_outlier': []
    }
    
    # We'll track metrics for ALL sequences, including removed ones
    metrics = {}
    
    # Function to update consensus and frequencies
    def update_consensus(records: List[SeqRecord]) -> Tuple[str, List[Dict[str, float]]]:
        if not records:
            return '', []
        temp_alignment = MultipleSeqAlignment(records)
        consensus, freqs = generate_consensus_sequence(temp_alignment, consensus_threshold)
        return str(consensus), freqs
    
    # Initialize consensus and frequencies - used only for initial metrics
    initial_consensus_seq, initial_frequencies = update_consensus(current_records)

    # Calculate initial metrics for all sequences
    for record in current_records:
        metrics[record.id] = analyze_sequence_metrics(
            record, initial_consensus_seq, initial_frequencies
    )
    
    # 1. Human COX1 similarity filtering
    if enable_human and current_records:
        kept, removed = filter_human_sequences(current_records, human_threshold)
        removed_records['human_similar'].extend(removed)
        current_records = kept
        
        # Always recalculate consensus after human filtering, before AT content comparison
        consensus_seq, frequencies = update_consensus(current_records)
        # Update metrics for remaining sequences
        for record in current_records:
            metrics[record.id].update(
                analyze_sequence_metrics(record, consensus_seq, frequencies)
            )
    
    # 2. AT content filtering with mode support
    if enable_at and current_records:
        kept, removed = filter_at_content(
            current_records,
            consensus_seq,
            threshold=at_threshold,
            mode=at_mode
        )
        removed_records['at_difference'].extend(removed)
        current_records = kept
        
        if removed:
            consensus_seq, frequencies = update_consensus(current_records)
            for record in current_records:
                metrics[record.id].update(
                    analyze_sequence_metrics(record, consensus_seq, frequencies)
                )
    
    # 3. Statistical outlier filtering
    if enable_outliers and current_records:
        kept, removed = filter_statistical_outliers(
            current_records, consensus_seq, frequencies, outlier_percentile
        )
        removed_records['statistical_outlier'].extend(removed)
        current_records = kept
        
        if removed:
            consensus_seq, frequencies = update_consensus(current_records)
            for record in current_records:
                metrics[record.id].update(
                    analyze_sequence_metrics(record, consensus_seq, frequencies)
                )
    
    # 4. Reference sequence filtering (if provided)
    if enable_reference and reference_seq and current_records:
        kept, removed = filter_reference_outliers(
            current_records, reference_seq, frequencies, outlier_percentile
        )
        removed_records['reference_outlier'].extend(removed)
        current_records = kept
        
        # Note: We don't recalculate consensus after reference filtering
        # as it's our last step and not needed for further filtering
    
    # Calculate reference metrics if reference sequence provided
    reference_metrics = None
    if reference_seq and enable_reference:
        reference_metrics = {
            'mean_deviation': np.mean([
                calculate_unweighted_deviation(str(record.seq), reference_seq)
                for record in current_records
            ]),
            'weighted_deviation': np.mean([
                calculate_weighted_deviation(str(record.seq), reference_seq, frequencies)
                for record in current_records
            ])
        }
    
    return FilterResult(
        kept_records=current_records,
        removed_records=removed_records,
        metrics=metrics,
        consensus_seq=consensus_seq,
        frequencies=frequencies,
        reference_metrics=reference_metrics
    )

def get_base_filename(filepath: str) -> str:
    """Extract base filename without extension."""
    basename = os.path.basename(filepath)
    for ext in ['.fasta', '.fas', '.fa']:
        if basename.lower().endswith(ext):
            return basename[:-len(ext)]
    return basename

def write_sequences_to_fasta(records: List[SeqRecord], filename: str) -> None:
    """Write sequences to FASTA file without line wrapping."""
    with open(filename, 'w') as handle:
        for record in records:
            handle.write(f">{record.id}\n{str(record.seq)}\n")

def write_metrics_report(filter_result: FilterResult, output_path: str) -> None:
    """
    Write comprehensive sequence metrics report including all filtering criteria.
    
    Args:
        filter_result: FilterResult containing sequence metrics and filtering results
        output_path: Path to write CSV report
    """
    with open(output_path, 'w', newline='') as csvfile:
        # Gather all possible metric fields
        metric_fields = set()
        for seq_metrics in filter_result.metrics.values():
            metric_fields.update(seq_metrics.keys())
        
        # Add standard fields
        base_fields = [
            'sequence_id',
            'length',
            'at_content',
            'human_similarity',
            'consensus_at',
            'at_difference',
            'unweighted_deviation',
            'weighted_deviation'
        ]
        
        # Add reference fields if available
        reference_fields = []
        if filter_result.reference_metrics:
            reference_fields = [
                'reference_unweighted_deviation',
                'reference_weighted_deviation',
                'reference_at_difference'
            ]
        
        # Position-specific metrics if available
        position_fields = []
        if filter_result.frequencies:
            position_fields = [
                'conservation_score',
                'gap_frequency'
            ]
        
        # Final fields list
        fields = base_fields + reference_fields + position_fields + ['removal_reason']
        
        writer = csv.writer(csvfile)
        writer.writerow(fields)
        
        # Create removal reason mapping
        removal_mapping = {}
        for reason, records in filter_result.removed_records.items():
            for record in records:
                removal_mapping[record.id] = reason
        
        # Write data for each sequence
        for seq_id, seq_metrics in filter_result.metrics.items():
            row = [seq_id]
            
            # Add base metrics
            for field in base_fields[1:]:  # Skip sequence_id
                value = seq_metrics.get(field, '')
                row.append(f"{value:.4f}" if isinstance(value, float) else str(value))
            
            # Add reference metrics
            if reference_fields:
                for field in reference_fields:
                    value = seq_metrics.get(field, '')
                    row.append(f"{value:.4f}" if isinstance(value, float) else str(value))
            
            # Add position-specific metrics
            if position_fields:
                for field in position_fields:
                    value = seq_metrics.get(field, '')
                    row.append(f"{value:.4f}" if isinstance(value, float) else str(value))
            
            # Add removal reason
            row.append(removal_mapping.get(seq_id, 'kept'))
            
            writer.writerow(row)

def process_fasta_file(input_file: str, 
                      output_dir: str,
                      reference_file: Optional[str] = None,
                      consensus_threshold: float = 0.5,
                      human_threshold: float = 0.95,
                      at_threshold: float = 0.1,
                      at_mode: str = 'absolute',
                      outlier_percentile: float = 90.0,
                      enable_human: bool = True,
                      enable_at: bool = True,
                      enable_outliers: bool = True,
                      enable_reference: bool = True) -> Dict[str, Any]:
    """
    Process a single FASTA file with all enabled filtering methods.
    
    Args:
        input_file: Path to input FASTA file
        output_dir: Directory for output files
        reference_file: Optional path to reference sequence file
        consensus_threshold: Threshold for consensus generation
        human_threshold: Human similarity threshold
        at_threshold: AT content difference threshold
        at_mode: AT content filtering mode ('absolute', 'higher', or 'lower')
        outlier_percentile: Percentile for outlier detection
        enable_human: Enable human similarity filtering
        enable_at: Enable AT content filtering
        enable_outliers: Enable statistical outlier detection
        enable_reference: Enable reference sequence filtering
    
    Returns:
        Dictionary containing processing statistics and filter results
    """
    # Set up output paths
    base_name = get_base_filename(input_file)
    output_paths = {
        'cleaned': os.path.join(output_dir, f"{base_name}_cleaned.fasta"),
        'removed': os.path.join(output_dir, f"{base_name}_removed_all.fasta"),
        'consensus': os.path.join(output_dir, f"{base_name}_consensus.fasta"),
        'metrics': os.path.join(output_dir, f"{base_name}_metrics.csv"),
        'log': os.path.join(output_dir, f"{base_name}_log.txt"),
        'ordered_annotated': os.path.join(output_dir, f"{base_name}_ordered_annotated.fasta")
    }
    
    # Create category-specific files for removed sequences
    removed_paths = {
        'human_similar': os.path.join(output_dir, f"{base_name}_removed_human.fasta"),
        'at_difference': os.path.join(output_dir, f"{base_name}_removed_at.fasta"),
        'statistical_outlier': os.path.join(output_dir, f"{base_name}_removed_outlier.fasta"),
        'reference_outlier': os.path.join(output_dir, f"{base_name}_removed_reference.fasta")
    }
    
    # Initialize statistics
    stats = {
        'input_sequences': 0,
        'kept_sequences': 0,
        'removed_sequences': defaultdict(int),
        'processing_time': None,
        'filter_result': None
    }
    
    start_time = datetime.now()
    
    with open(output_paths['log'], 'w') as log_file:
        try:
            # Log start time and parameters
            log_file.write(f"Processing started at: {start_time}\n")
            log_file.write(f"Input file: {input_file}\n\n")
            log_file.write("Parameters:\n")
            log_file.write(f"- Consensus threshold: {consensus_threshold}\n")
            log_file.write(f"- Human threshold: {human_threshold}\n")
            log_file.write(f"- AT threshold: {at_threshold}\n")
            log_file.write(f"- AT mode: {at_mode}\n")
            log_file.write(f"- Outlier percentile: {outlier_percentile}\n")
            log_file.write("Enabled filters:\n")
            log_file.write(f"- Human similarity: {'Yes' if enable_human else 'No'}\n")
            log_file.write(f"- AT content: {'Yes' if enable_at else 'No'}\n")
            log_file.write(f"- Statistical outliers: {'Yes' if enable_outliers else 'No'}\n")
            log_file.write(f"- Reference comparison: {'Yes' if enable_reference else 'No'}\n\n")
            
            # Read input alignment
            alignment = AlignIO.read(input_file, "fasta")
            stats['input_sequences'] = len(alignment)
            log_file.write(f"Input sequences: {stats['input_sequences']}\n")
            
            # Get reference sequence if available
            reference_seq = None
            if reference_file and enable_reference:
                try:
                    reference_seq = get_reference_sequence(reference_file)
                    log_file.write(f"Using reference file: {reference_file}\n")
                except Exception as e:
                    log_file.write(f"Warning: Could not read reference file - {str(e)}\n")
                    enable_reference = False
            
            # Apply filters
            filter_result = apply_filters(
                alignment,
                consensus_threshold=consensus_threshold,
                human_threshold=human_threshold,
                at_threshold=at_threshold,
                at_mode=at_mode,
                outlier_percentile=outlier_percentile,
                reference_seq=reference_seq,
                enable_human=enable_human,
                enable_at=enable_at,
                enable_outliers=enable_outliers,
                enable_reference=enable_reference
            )
            
            stats['filter_result'] = filter_result
            
            # Write outputs
            if filter_result.kept_records:
                # Sort kept sequences by start position before writing
                sorted_kept_records = sort_sequences_by_position(filter_result.kept_records)
                write_sequences_to_fasta(sorted_kept_records, output_paths['cleaned'])
                log_file.write(f"Wrote {len(sorted_kept_records)} kept sequences to: {output_paths['cleaned']}\n")
               
                # Write consensus sequence
                consensus_record = SeqRecord(
                    Seq(filter_result.consensus_seq),
                    id=f"{base_name}_consensus",
                    description=f"c{consensus_threshold}_h{human_threshold}_a{at_threshold}_p{outlier_percentile}"
                )
                write_sequences_to_fasta([consensus_record], output_paths['consensus'])
                log_file.write(f"Wrote consensus sequence to: {output_paths['consensus']}\n")
            
            # Write removed sequences by category
            all_removed = []
            for reason, records in filter_result.removed_records.items():
                if records:
                    write_sequences_to_fasta(records, removed_paths[reason])
                    all_removed.extend(records)
                    stats['removed_sequences'][reason] = len(records)
                    log_file.write(f"Wrote {len(records)} {reason} sequences to: {removed_paths[reason]}\n")
            
            # Write combined removed sequences
            if all_removed:
                write_sequences_to_fasta(all_removed, output_paths['removed'])
                log_file.write(f"Wrote {len(all_removed)} total removed sequences to: {output_paths['removed']}\n")
            
             # Write new ordered and annotated alignment file
            write_ordered_annotated_alignment(
                filter_result.kept_records,
                filter_result.removed_records,
                output_paths['ordered_annotated']
            )
            log_file.write(f"Wrote ordered and annotated alignment to: {output_paths['ordered_annotated']}\n")
            
            # Write metrics report
            write_metrics_report(filter_result, output_paths['metrics'])
            log_file.write(f"Wrote metrics report to: {output_paths['metrics']}\n")
            
            # Update statistics
            stats['kept_sequences'] = len(filter_result.kept_records)
            stats['processing_time'] = (datetime.now() - start_time).total_seconds()
            
            # Log final statistics
            log_file.write("\nProcessing Results:\n")
            log_file.write(f"Input sequences: {stats['input_sequences']}\n")
            log_file.write(f"Kept sequences: {stats['kept_sequences']}\n")
            for reason, count in stats['removed_sequences'].items():
                log_file.write(f"Removed ({reason}): {count}\n")
            log_file.write(f"\nProcessing completed in {stats['processing_time']:.2f} seconds\n")
            
        except Exception as e:
            error_msg = f"Error processing {input_file}: {str(e)}"
            log_file.write(f"\nERROR: {error_msg}\n")
            log_file.write(f"Traceback:\n{traceback.format_exc()}\n")
            raise RuntimeError(error_msg)
    
    return stats

def process_directory(args: argparse.Namespace) -> None:
    """
    Process all FASTA files in the input directory.
    
    Args:
        args: Parsed command line arguments
    """
    total_stats = {
        'processed_files': 0,
        'total_sequences': 0,
        'kept_sequences': 0,
        'removed_sequences': {
            'human_similar': 0,
            'at_difference': 0,
            'statistical_outlier': 0
        }
    }
    
    # Process each FASTA file
    for filename in os.listdir(args.input_dir):
        if filename.lower().endswith(('.fasta', '.fas', '.fa')):
            input_path = os.path.join(args.input_dir, filename)
            
            # Get reference file path if applicable
            reference_file = None
            if args.reference_dir:
                base_name = get_base_filename(filename)
                reference_path = os.path.join(args.reference_dir, f"{base_name}_reference.fasta")
                if os.path.exists(reference_path):
                    reference_file = reference_path
            
            try:
                stats = process_fasta_file(
                    input_path,
                    args.output_dir,
                    reference_file=reference_file,
                    consensus_threshold=args.consensus_threshold,
                    human_threshold=args.human_threshold,
                    at_threshold=args.at_difference,
                    at_mode=args.at_mode,
                    outlier_percentile=args.percentile_threshold,
                    enable_human=not args.disable_human,
                    enable_at=not args.disable_at,
                    enable_outliers=not args.disable_outliers
                )
                
                # Update total statistics
                total_stats['processed_files'] += 1
                total_stats['total_sequences'] += stats['input_sequences']
                total_stats['kept_sequences'] += stats['kept_sequences']
                for reason, count in stats['removed_sequences'].items():
                    total_stats['removed_sequences'][reason] += count
                
                print(f"Processed {filename}")
                print(f"- Input sequences: {stats['input_sequences']}")
                print(f"- Kept sequences: {stats['kept_sequences']}")
                for reason, count in stats['removed_sequences'].items():
                    print(f"- Removed ({reason}): {count}")
                print()
                
            except Exception as e:
                print(f"Error processing {filename}: {str(e)}")
                continue
    
    # Print final summary
    print("\nProcessing Summary:")
    print(f"Processed {total_stats['processed_files']} files")
    print(f"Total sequences: {total_stats['total_sequences']}")
    print(f"Total kept: {total_stats['kept_sequences']}")
    for reason, count in total_stats['removed_sequences'].items():
        print(f"Total removed ({reason}): {count}")

if __name__ == "__main__":
    try:
        # Parse command line arguments
        args = parse_arguments()
        
        # Process all files in the input directory
        process_directory(args)
        
    except Exception as e:
        print(f"Error: {str(e)}")
        sys.exit(1)
