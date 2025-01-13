import os
import csv
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import substitution_matrices, PairwiseAligner
from Bio.Align import MultipleSeqAlignment
from Bio.motifs import Motif
from datetime import datetime
import tempfile
import traceback
from typing import List, Tuple, Dict, Optional

def parse_arguments():
    """Set up command line argument parsing with clear help text and validation."""
    parser = argparse.ArgumentParser(
        description='Process FASTA files to filter sequences based on target vs contaminant matching.',
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
        '-n', '--nontarget',
        required=True,
        help='Path to FASTA file containing non-target contaminant sequences'
    )
    
    parser.add_argument(
        '-t', '--targets_dir',
        required=True,
        help='Directory containing target sequence FASTA files'
    )
    
    parser.add_argument(
        '-c', '--consensus_threshold',
        type=float,
        default=0.7,
        help='Consensus threshold (0.0-1.0). Minimum frequency required for a base to be included in consensus'
    )
    
    parser.add_argument(
        '--translation_table',
        type=int,
        default=1,
        help='Genetic code table to use for translation (default: 1 Standard). Use 5 for invertebrate mitochondrial code'
    )
    
    args = parser.parse_args()
    
    if not 0 <= args.consensus_threshold <= 1:
        parser.error("Consensus threshold must be between 0.0 and 1.0")
    
    # Validate translation table
    valid_tables = [1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 16, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 33]
    if args.translation_table not in valid_tables:
        parser.error(f"Translation table {args.translation_table} is not valid. Valid tables are: {valid_tables}")
    
    # Validate directories and files exist
    if not os.path.isdir(args.input_dir):
        parser.error(f"Input directory does not exist: {args.input_dir}")
    
    if not os.path.isdir(args.targets_dir):
        parser.error(f"Targets directory does not exist: {args.targets_dir}")
        
    if not os.path.isfile(args.nontarget):
        parser.error(f"Contaminants file does not exist: {args.nontarget}")
    
    return args

def create_protein_aligner() -> PairwiseAligner:
    """Create a protein aligner with BLOSUM62 settings."""
    aligner = PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    aligner.mode = 'global'
    return aligner

def remove_gaps(sequence: str) -> str:
    """Remove gap characters from sequence."""
    return str(sequence).replace('-', '')

def translate_frames(dna_sequence: str, translation_table: int = 1) -> List[str]:
    """
    Translate DNA sequence in all 3 reading frames using specified translation table.
    
    Args:
        dna_sequence: DNA sequence to translate
        translation_table: Genetic code table number (default: 1 Standard)
            Common tables:
            1 = Standard
            2 = Vertebrate Mitochondrial
            5 = Invertebrate Mitochondrial
            11 = Bacterial and Plant Plastid
    """
    seq = Seq(dna_sequence)
    frames = []
    
    # Forward frames
    for i in range(3):
        # Remove any trailing incomplete codons
        truncated_len = len(dna_sequence[i:]) - len(dna_sequence[i:]) % 3
        frame = seq[i:i + truncated_len].translate(table=translation_table)
        frames.append(str(frame))
    
    return frames

def calculate_alignment_stats(alignment) -> Dict[str, float]:
    """Calculate identity, similarity, and gap percentages from alignment."""
    seq1, seq2 = str(alignment[0]), str(alignment[1])
    length = len(seq1)
    matches = sum(1 for i in range(length) if seq1[i] == seq2[i])
    gaps = sum(1 for c in seq1 if c == '-') + sum(1 for c in seq2 if c == '-')
    
    # Calculate similarity using BLOSUM62 matrix
    matrix = substitution_matrices.load("BLOSUM62")
    similar = 0
    for i in range(length):
        if seq1[i] != '-' and seq2[i] != '-':
            if matrix.get((seq1[i], seq2[i]), -4) > 0:
                similar += 1
    
    return {
        'identity': (matches / length) * 100,
        'similarity': (similar / length) * 100,
        'gaps': (gaps / (2 * length)) * 100
    }

def format_alignment_output(alignment, query_name: str, target_name: str) -> str:
    """Format alignment for logging with match/mismatch markers."""
    aligned_seq1, aligned_seq2 = str(alignment[0]), str(alignment[1])
    
    # Create match line
    match_line = ''
    for i in range(len(aligned_seq1)):
        if aligned_seq1[i] == aligned_seq2[i]:
            match_line += '|'
        elif aligned_seq1[i] == '-' or aligned_seq2[i] == '-':
            match_line += ' '
        else:
            # Check for similar amino acids using BLOSUM62
            matrix = substitution_matrices.load("BLOSUM62")
            if matrix.get((aligned_seq1[i], aligned_seq2[i]), -4) > 0:
                match_line += ':'
            else:
                match_line += ' '
    
    # Calculate statistics
    stats = calculate_alignment_stats(alignment)
    
    # Format output
    output = [
        f"{target_name:<15} {aligned_seq1}",
        f"{'':15} {match_line}",
        f"{query_name:<15} {aligned_seq2}",
        f"Statistics:",
        f"  Identity:   {stats['identity']:.1f}%",
        f"  Similarity: {stats['similarity']:.1f}%",
        f"  Gaps:       {stats['gaps']:.1f}%"
    ]
    
    return '\n'.join(output)

def align_sequences(query_seq: str, reference_seq: str, aligner: PairwiseAligner) -> Tuple[float, str]:
    """Perform pairwise alignment between query and reference sequences."""
    alignments = aligner.align(query_seq, reference_seq)
    alignment = alignments[0]  # Get first (best) alignment
    return alignment.score, alignment

def get_best_frame_match(frames: List[str], reference_seq: str, query_name: str, 
                        reference_name: str, aligner: PairwiseAligner, log_file) -> Tuple[str, float, str]:
    """Find the best matching frame and its alignment score."""
    best_score = float('-inf')
    best_frame = None
    best_alignment_output = ""
    
    with open(log_file, 'a') as f:
        f.write(f"\nAligning frames against {reference_name}:\n")
        f.write("=" * 60 + "\n")
        
        for i, frame in enumerate(frames, 1):
            score, alignment = align_sequences(frame, reference_seq, aligner)
            alignment_output = format_alignment_output(alignment, f"Frame {i}", reference_name)
            
            f.write(f"\nFrame {i} (Score: {score:.2f}):\n")
            f.write(alignment_output + "\n")
            
            if score > best_score:
                best_score = score
                best_frame = frame
                best_alignment_output = alignment_output
        
        f.write(f"\nSelected Frame {frames.index(best_frame) + 1} as best match (Score: {best_score:.2f})\n")
    
    return best_frame, best_score, best_alignment_output

def get_best_contaminant_match(frames: List[str], contaminant_seqs: List[SeqRecord], 
                              query_name: str, aligner: PairwiseAligner, log_file) -> Tuple[str, float]:
    """Find the best matching contaminant and its alignment score."""
    best_score = float('-inf')
    best_contaminant = None
    
    for contaminant in contaminant_seqs:
        frame, score, alignment_output = get_best_frame_match(
            frames, str(contaminant.seq), query_name, contaminant.id, aligner, log_file
        )
        if score > best_score:
            best_score = score
            best_contaminant = contaminant.id
    
    return best_contaminant, best_score

def process_single_read(read: SeqRecord, target_seq: SeqRecord, 
                       contaminant_seqs: List[SeqRecord], log_file: str,
                       translation_table: int = 1) -> Dict:
    """Process a single read with detailed logging."""
    aligner = create_protein_aligner()
    
    with open(log_file, 'a') as f:
        f.write(f"\n{'='*80}\n")
        f.write(f"Processing Read: {read.id}\n")
        f.write(f"{'='*80}\n\n")
        
        # Add translation table info to log
        f.write(f"Using translation table: {translation_table}\n\n")
        
        # Log original sequence
        f.write("Original DNA Sequence (with gaps):\n")
        f.write(str(read.seq) + "\n\n")
        
        # Remove gaps and translate
        ungapped = remove_gaps(read.seq)
        f.write("Ungapped DNA Sequence:\n")
        f.write(ungapped + "\n\n")
        
        # Translate frames
        frames = translate_frames(ungapped, translation_table)
        f.write("=== Frame Translations ===\n")
        for i, frame in enumerate(frames, 1):
            f.write(f"Frame {i}:\n{frame}\n\n")
        
        # Process against target
        f.write("=== Target Sequence Comparison ===\n")
        best_target_frame, target_score, target_alignment = get_best_frame_match(
            frames, str(target_seq.seq), read.id, "Target", aligner, log_file
        )
        
        # Process against contaminants
        f.write("\n=== Contaminant Sequence Comparisons ===\n")
        best_contaminant, contaminant_score = get_best_contaminant_match(
            frames, contaminant_seqs, read.id, aligner, log_file
        )
        
        # Log final decision
        f.write("\n=== Final Decision ===\n")
        decision = target_score > contaminant_score
        f.write(f"Best Target Score: {target_score:.2f}\n")
        f.write(f"Best Contaminant Score: {contaminant_score:.2f} (Contaminant: {best_contaminant})\n")
        f.write(f"Decision: {'KEPT' if decision else 'REMOVED'}\n")
        f.write(f"{'='*80}\n\n")
    
    return {
        'keep': decision,
        'target_score': target_score,
        'contaminant_score': contaminant_score,
        'best_contaminant': best_contaminant
    }

def write_sequence_report(sequences, results, output_path):
    """Write sequence filtering results to CSV file."""
    with open(output_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow([
            'Sequence Name',
            'Status',
            'Target Score',
            'Contaminant Score',
            'Best Matching Contaminant'
        ])
        
        for record, result in zip(sequences, results):
            writer.writerow([
                record.id,
                'kept' if result['keep'] else 'removed',
                f"{result['target_score']:.3f}",
                f"{result['contaminant_score']:.3f}",
                result['best_contaminant']
            ])

def generate_consensus(records, identity_threshold=0.7, log_file=None):
    """Generate consensus sequence from aligned records with detailed logging."""
    if log_file:
        with open(log_file, 'a') as f:
            f.write("\n=== Consensus Generation ===\n")
            f.write(f"Starting consensus generation at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Number of input sequences: {len(records)}\n")
            f.write(f"Identity threshold: {identity_threshold}\n\n")
    
    try:
        if not records:
            if log_file:
                with open(log_file, 'a') as f:
                    f.write("Error: No sequences provided for consensus generation\n")
            return None
        
        # Get sequence length from first record
        seq_length = len(records[0].seq)
        consensus_seq = []
        
        if log_file:
            with open(log_file, 'a') as f:
                f.write(f"Sequence length: {seq_length} bases\n")
                f.write("Beginning position-by-position analysis...\n\n")
        
        # Process each position in the alignment
        for i in range(seq_length):
            try:
                # Get all bases at this position
                column = [record.seq[i] for record in records]
                
                # Track position details for logging
                position_details = {
                    'position': i + 1,
                    'total_bases': len(column),
                    'gaps': column.count('-'),
                    'base_counts': {}
                }
                
                # If this position is a gap in all sequences, keep it as a gap
                if all(base == '-' for base in column):
                    consensus_seq.append('-')
                    if log_file:
                        with open(log_file, 'a') as f:
                            f.write(f"Position {i+1}: All gaps - keeping gap in consensus\n")
                    continue
                
                # Remove gaps for calculating base frequencies
                bases = [base for base in column if base != '-']
                
                if not bases:
                    consensus_seq.append('N')
                    if log_file:
                        with open(log_file, 'a') as f:
                            f.write(f"Position {i+1}: No non-gap bases found - using 'N'\n")
                    continue
                
                # Count frequencies of bases
                base_counts = {}
                for base in bases:
                    base_counts[base] = base_counts.get(base, 0) + 1
                position_details['base_counts'] = base_counts
                
                # Find most common base and its frequency
                most_common_base = max(base_counts.items(), key=lambda x: x[1])
                frequency = most_common_base[1] / len(bases)
                
                # Add to consensus if frequency meets threshold
                if frequency >= identity_threshold:
                    consensus_seq.append(most_common_base[0])
                    if log_file:
                        with open(log_file, 'a') as f:
                            f.write(f"Position {i+1}: Selected {most_common_base[0]} "
                                  f"(frequency: {frequency:.2f}, counts: {base_counts})\n")
                else:
                    consensus_seq.append('N')
                    if log_file:
                        with open(log_file, 'a') as f:
                            f.write(f"Position {i+1}: Below threshold - using 'N' "
                                  f"(frequency: {frequency:.2f}, counts: {base_counts})\n")
            
            except Exception as e:
                # Handle errors at individual positions
                consensus_seq.append('N')
                if log_file:
                    with open(log_file, 'a') as f:
                        f.write(f"Error at position {i+1}: {str(e)} - using 'N'\n")
        
        # Create consensus sequence record
        consensus_sequence = ''.join(consensus_seq)
        consensus_record = SeqRecord(
            Seq(consensus_sequence),
            id="consensus",
            description=f"consensus_identity_{identity_threshold}"
        )
        
        # Log completion and statistics
        if log_file:
            with open(log_file, 'a') as f:
                f.write("\nConsensus Generation Summary:\n")
                f.write("=" * 30 + "\n")
                f.write(f"Total length: {len(consensus_sequence)}\n")
                f.write(f"Number of gaps: {consensus_sequence.count('-')}\n")
                f.write(f"Number of ambiguous positions (N): {consensus_sequence.count('N')}\n")
                f.write(f"Completed at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        
        return consensus_record
    
    except Exception as e:
        # Handle any unexpected errors during consensus generation
        if log_file:
            with open(log_file, 'a') as f:
                f.write(f"\nError during consensus generation: {str(e)}\n")
                f.write(f"Traceback:\n{traceback.format_exc()}\n")
        raise

def process_directory(input_dir: str, output_dir: str, targets_dir: str, 
                     contaminants_path: str, consensus_threshold: float = 0.7,
                     translation_table: int = 1):
    """Process all FASTA files in a directory."""
    os.makedirs(output_dir, exist_ok=True)
    
    # Load contaminant sequences once
    contaminant_seqs = list(SeqIO.parse(contaminants_path, "fasta"))
    processed_count = 0
    
    for filename in os.listdir(input_dir):
        if filename.lower().endswith(('.fasta', '.fas', '.fa')):
            input_path = os.path.join(input_dir, filename)
            base_name = os.path.splitext(filename)[0]
            
            try:
                # Find corresponding target file
                target_path = None
                for ext in ['.fasta', '.fas', '.fa']:
                    possible_path = os.path.join(targets_dir, f"{base_name.replace('_align', '')}{ext}")
                    if os.path.exists(possible_path):
                        target_path = possible_path
                        break
                
                if not target_path:
                    raise FileNotFoundError(f"No target file found for {filename}")
                
                # Create log file
                log_path = os.path.join(output_dir, f"{base_name}_detailed.log")
                with open(log_path, 'w') as f:
                    f.write(f"Detailed processing log for {filename}\n")
                    f.write(f"Generated at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
                    f.write(f"Using translation table: {translation_table}\n")
                    f.write("="*80 + "\n\n")
                
                # Process sequences
                try:
                    target = next(SeqIO.parse(target_path, "fasta"))
                except StopIteration:
                    print(f"Error: Target file {target_path} is empty")
                    continue
                except Exception as e:
                    print(f"Error reading target file {target_path}: {str(e)}")
                    continue
                
                try:
                    alignment = list(SeqIO.parse(input_path, "fasta"))
                    if not alignment:
                        print(f"Error: Input file {input_path} is empty")
                        continue
                except Exception as e:
                    print(f"Error reading input file {input_path}: {str(e)}")
                    continue
                
                kept_records = []
                removed_records = []
                results = []
                
                for record in alignment:
                    result = process_single_read(record, target, contaminant_seqs, 
                                              log_path, translation_table)
                    results.append(result)
                    
                    if result['keep']:
                        kept_records.append(record)
                    else:
                        removed_records.append(record)

                # Generate consensus from kept records
                consensus_record = None
                try:
                    if kept_records:
                        consensus_record = generate_consensus(
                            kept_records, 
                            consensus_threshold,
                            log_path
                        )
                    else:
                        with open(log_path, 'a') as f:
                            f.write("\nSkipping consensus generation: No sequences kept after filtering\n")
                except Exception as e:
                    print(f"Warning: Failed to generate consensus for {filename}: {str(e)}")
                    with open(log_path, 'a') as f:
                        f.write(f"\nError during consensus generation:\n{str(e)}\n")
                        f.write(f"Traceback:\n{traceback.format_exc()}\n")
                
                # Write outputs with error handling
                try:
                    cleaned_path = os.path.join(output_dir, f"{base_name}_cleaned.fasta")
                    removed_path = os.path.join(output_dir, f"{base_name}_removed.fasta")
                    consensus_path = os.path.join(output_dir, f"{base_name}_consensus.fasta")
                    report_path = os.path.join(output_dir, f"{base_name}_report.csv")
                    
                    write_fasta_single_line(kept_records, cleaned_path)
                    write_fasta_single_line(removed_records, removed_path)
                    if consensus_record:
                        write_fasta_single_line([consensus_record], consensus_path)
                    write_sequence_report(alignment, results, report_path)
                except Exception as e:
                    print(f"Error writing output files for {filename}: {str(e)}")
                    with open(log_path, 'a') as f:
                        f.write(f"\nError writing output files:\n{str(e)}\n")
                        f.write(f"Traceback:\n{traceback.format_exc()}\n")
                    continue
                
                processed_count += 1
                
            except FileNotFoundError as e:
                print(f"Error processing {filename}: {str(e)}")
                continue
            except Exception as e:
                print(f"Unexpected error processing {filename}: {str(e)}")
                continue
    
    if processed_count == 0:
        print(f"No FASTA files processed in {input_dir}")
    else:
        print(f"Finished processing {processed_count} FASTA files")

def write_fasta_single_line(records: List[SeqRecord], filename: str):
    """Write sequences to FASTA file without line wrapping."""
    with open(filename, 'w') as handle:
        for record in records:
            handle.write(f">{record.id} {record.description}\n")
            handle.write(f"{str(record.seq)}\n")

if __name__ == "__main__":
    # Set up argument parsing
    args = parse_arguments()
    
    # Process all files in the input directory
    try:
        process_directory(
            args.input_dir,
            args.output_dir,
            args.targets_dir,
            args.nontarget,
            args.consensus_threshold,
            args.translation_table
        )
    except Exception as e:
        print(f"Error during execution: {str(e)}")
        raise
