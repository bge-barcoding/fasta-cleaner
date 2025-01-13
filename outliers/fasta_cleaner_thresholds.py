import os
import csv
import argparse
from Bio import SeqIO
from collections import Counter
from typing import List, Dict, Tuple
import numpy as np
from datetime import datetime

def parse_arguments():
    """Set up command line argument parsing with clear help text."""
    parser = argparse.ArgumentParser(
        description='Analyze FASTA alignment files to identify outlier sequences.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        '-i', '--input_dir',
        required=True,
        help='Directory containing input alignment FASTA files'
    )
    
    parser.add_argument(
        '-o', '--output_dir',
        required=True,
        help='Directory where output files will be saved'
    )
    
    parser.add_argument(
        '-r', '--reference_dir',
        help='Optional directory containing reference FASTA files (named same as input files but with "_reference.fasta" suffix)'
    )

    parser.add_argument(
        '-c', '--consensus_threshold',
        type=float,
        default=0.7,
        help='Threshold for consensus sequence generation (between 0 and 1)'
    )
    
    args = parser.parse_args()
    
    # Validate directories exist
    if not os.path.isdir(args.input_dir):
        parser.error(f"Input directory does not exist: {args.input_dir}")
    
    if args.reference_dir and not os.path.isdir(args.reference_dir):
        parser.error(f"Reference directory does not exist: {args.reference_dir}")
    
    # Validate consensus threshold
    if not 0 < args.consensus_threshold <= 1:
        parser.error("Consensus threshold must be between 0 and 1")
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)
    
    return args

def get_reference_sequence(reference_file: str) -> str:
    """Extract reference sequence from a FASTA file."""
    try:
        references = list(SeqIO.parse(reference_file, "fasta"))
        if not references:
            raise ValueError("No sequences found in reference file")
        if len(references) > 1:
            print(f"Warning: Multiple sequences found in {reference_file}, using first one")
        return str(references[0].seq)
    except Exception as e:
        raise ValueError(f"Error reading reference file: {str(e)}")

def calculate_consensus_and_frequencies(sequences: List[str]) -> Tuple[str, List[Dict[str, float]]]:
    """Calculate consensus sequence and position-specific frequencies."""
    if not sequences:
        raise ValueError("No sequences provided")
    
    align_length = len(sequences[0])
    if not all(len(seq) == align_length for seq in sequences):
        raise ValueError("Sequences have different lengths - alignment required")
    
    frequencies = calculate_position_frequencies(sequences)
    consensus = ''
    for pos_freqs in frequencies:
        if pos_freqs:
            consensus += max(pos_freqs.items(), key=lambda x: x[1])[0]
        else:
            consensus += '-'
    
    return consensus, frequencies

def calculate_weighted_deviation_score(sequence: str, reference: str, 
                                    frequencies: List[Dict[str, float]]) -> float:
    """Calculate weighted deviation score between sequence and reference."""
    if len(sequence) != len(reference) != len(frequencies):
        raise ValueError("Sequence, reference, and frequencies must all have same length")
    
    total_score = 0.0
    total_weight = 0.0
    
    for seq_res, ref_res, pos_freqs in zip(sequence, reference, frequencies):
        if seq_res != '-' and ref_res != '-':
            conservation_weight = pos_freqs.get(ref_res, 0)
            total_weight += conservation_weight
            
            if seq_res != ref_res:
                total_score += conservation_weight
    
    return total_score / total_weight if total_weight > 0 else 0.0

def calculate_unweighted_deviation_score(sequence: str, reference: str) -> float:
    """Calculate unweighted deviation score between sequence and reference."""
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

def calculate_position_frequencies(sequences: List[str]) -> List[Dict[str, float]]:
    """Calculate residue frequencies at each position in the alignment."""
    if not sequences:
        raise ValueError("No sequences provided")
    
    align_length = len(sequences[0])
    position_counts = []
    
    for i in range(align_length):
        residues = [seq[i] for seq in sequences if seq[i] != '-']
        if residues:
            counts = Counter(residues)
            total = sum(counts.values())
            frequencies = {res: count/total for res, count in counts.items()}
        else:
            frequencies = {}
        position_counts.append(frequencies)
    
    return position_counts

def generate_consensus_sequence(sequences: List[str], threshold: float) -> str:
    """Generate consensus sequence using specified threshold."""
    if not sequences:
        return ""
    
    align_length = len(sequences[0])
    consensus = []
    
    for i in range(align_length):
        # Get residues at this position
        residues = [seq[i] for seq in sequences if seq[i] != '-']
        if not residues:
            consensus.append('-')
            continue
            
        # Count frequencies
        counts = Counter(residues)
        total = len(residues)
        
        # Find most common residue and its frequency
        most_common = counts.most_common(1)[0]
        if most_common[1] / total >= threshold:
            consensus.append(most_common[0])
        else:
            consensus.append('-')
    
    return ''.join(consensus)

def analyze_alignment_file(alignment_file: str, reference_sequence: str = None) -> List[Tuple[str, float, float, float, float]]:
    """Analyze alignment file and calculate deviation scores."""
    sequences = []
    seq_ids = []
    for record in SeqIO.parse(alignment_file, "fasta"):
        sequences.append(str(record.seq))
        seq_ids.append(record.id)
    
    consensus, frequencies = calculate_consensus_and_frequencies(sequences)
    
    scores = []
    for seq_id, sequence in zip(seq_ids, sequences):
        consensus_unweighted = calculate_unweighted_deviation_score(sequence, consensus)
        consensus_weighted = calculate_weighted_deviation_score(sequence, consensus, frequencies)
        
        if reference_sequence:
            reference_unweighted = calculate_unweighted_deviation_score(sequence, reference_sequence)
            reference_weighted = calculate_weighted_deviation_score(
                sequence, 
                reference_sequence,
                frequencies
            )
        else:
            reference_unweighted = None
            reference_weighted = None
            
        scores.append((
            seq_id, 
            consensus_unweighted, 
            consensus_weighted,
            reference_unweighted,
            reference_weighted
        ))
    
    return scores, sequences, seq_ids

def write_sequences_to_fasta(sequences: List[str], seq_ids: List[str], output_file: str):
    """Write sequences to a FASTA file."""
    with open(output_file, 'w') as f:
        for seq_id, sequence in zip(seq_ids, sequences):
            f.write(f">{seq_id}\n{sequence}\n")

def process_directory(input_dir: str, output_dir: str, reference_dir: str = None, consensus_threshold: float = 0.7):
    """Process all FASTA files in the input directory."""
    for filename in os.listdir(input_dir):
        if filename.lower().endswith(('.fasta', '.fas', '.fa')):
            input_path = os.path.join(input_dir, filename)
            base_name = os.path.splitext(filename)[0]
            output_path = os.path.join(output_dir, f"{base_name}_scores.csv")
            cleaned_path = os.path.join(output_dir, f"{base_name}_cleaned.fasta")
            removed_path = os.path.join(output_dir, f"{base_name}_removed.fasta")
            log_path = os.path.join(output_dir, f"{base_name}_log.txt")
            
            with open(log_path, 'w') as log_file:
                log_file.write(f"Processing started at: {datetime.now()}\n")
                log_file.write(f"Input file: {input_path}\n")
                
                try:
                    # Get reference sequence if available
                    reference_sequence = None
                    if reference_dir:
                        reference_file = os.path.join(reference_dir, f"{base_name}_reference.fasta")
                        if os.path.exists(reference_file):
                            reference_sequence = get_reference_sequence(reference_file)
                            log_file.write(f"Using reference file: {reference_file}\n")
                        else:
                            log_file.write(f"No reference file found for {filename}\n")
                    
                    # Analyze alignment file
                    scores, sequences, seq_ids = analyze_alignment_file(input_path, reference_sequence)
                    
                    # Calculate thresholds using non-zero scores
                    non_zero_scores = [score[1] for score in scores if score[1] > 0]  # Using unweighted scores
                    mean_threshold = np.mean(non_zero_scores)
                    percentile_75_threshold = np.percentile(non_zero_scores, 75)
                    
                    log_file.write(f"\nScore thresholds:\n")
                    log_file.write(f"Mean threshold: {mean_threshold:.4f}\n")
                    log_file.write(f"75th percentile threshold: {percentile_75_threshold:.4f}\n")
                    
                    # Split sequences based on 75th percentile threshold
                    cleaned_sequences = []
                    cleaned_ids = []
                    removed_sequences = []
                    removed_ids = []
                    
                    for (seq_id, unweighted, weighted, ref_unweighted, ref_weighted), sequence in zip(scores, sequences):
                        if unweighted <= percentile_75_threshold:
                            cleaned_sequences.append(sequence)
                            cleaned_ids.append(seq_id)
                        else:
                            removed_sequences.append(sequence)
                            removed_ids.append(seq_id)
                    
                    # Write filtered sequences
                    write_sequences_to_fasta(cleaned_sequences, cleaned_ids, cleaned_path)
                    write_sequences_to_fasta(removed_sequences, removed_ids, removed_path)
                    
                    log_file.write(f"\nSequence filtering results:\n")
                    log_file.write(f"Kept sequences: {len(cleaned_sequences)}\n")
                    log_file.write(f"Removed sequences: {len(removed_sequences)}\n")
                    
                    # Generate consensus from cleaned sequences
                    consensus = generate_consensus_sequence(cleaned_sequences, consensus_threshold)
                    consensus_path = os.path.join(output_dir, f"{base_name}_consensus.fasta")
                    with open(consensus_path, 'w') as f:
                        f.write(f">{base_name}_consensus\n{consensus}\n")
                    
                    log_file.write(f"\nConsensus sequence generated with threshold {consensus_threshold}\n")
                    
                    # Prepare fate information for each sequence
                    sequence_fates = {}
                    for seq_id in seq_ids:
                        if seq_id in cleaned_ids:
                            sequence_fates[seq_id] = "kept"
                        else:
                            sequence_fates[seq_id] = "removed"

                    # Write scores to CSV with fate information
                    with open(output_path, 'w', newline='') as csvfile:
                        writer = csv.writer(csvfile)
                        headers = [
                            'Sequence ID',
                            'Consensus Unweighted Score',
                            'Consensus Weighted Score',
                            'Fate'  # New column
                        ]
                        if reference_sequence:
                            headers.extend([
                                'Reference Unweighted Score',
                                'Reference Weighted Score'
                            ])
                        writer.writerow(headers)
                        
                        for score_tuple in scores:
                            seq_id = score_tuple[0]
                            row = [seq_id]
                            row.extend(f"{score:.4f}" for score in score_tuple[1:3])
                            row.append(sequence_fates[seq_id])  # Add fate
                            if reference_sequence:
                                row.extend(f"{score:.4f}" for score in score_tuple[3:])
                            writer.writerow(row)
                    
                    log_file.write("\nProcessing completed successfully\n")
                    
                except Exception as e:
                    error_message = f"Error processing {filename}: {str(e)}"
                    print(error_message)
                    log_file.write(f"\nERROR: {error_message}\n")
                    continue

if __name__ == "__main__":
    args = parse_arguments()
    process_directory(args.input_dir, args.output_dir, args.reference_dir, args.consensus_threshold)
